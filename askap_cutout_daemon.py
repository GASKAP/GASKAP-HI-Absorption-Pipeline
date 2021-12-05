#!/usr/bin/env python -u

# Daemon process to manage the production of cutouts from an ASKAP scheduling block

# Author James Dempsey
# Date 27 Jan 2020

import argparse
import csv
import datetime
import glob

import os
import subprocess
import sys
import time
from pathlib import Path

from astropy.io.votable import parse_single_table


class CommandFailedError(Exception):
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


def parseargs():
    """
    Parse the command line arguments
    :return: An args map with the parsed arguments
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     description="Daemon process to manage the production of cutouts for each component from an ASKAP scheduing block")
    parser.add_argument("-d", "--delay", help="Number of seconds to pause between scans for completed jobs",
                        type=int, default=30)
    parser.add_argument("-s", "--sbid", help="The id of the ASKAP scheduling block to be processed",
                        type=int, required=True)
    parser.add_argument("--status_folder", help="The status folder which will contain the job completion or failed files",
                        default='status')
    parser.add_argument("-f", "--filename", help="The name of the csv format file listing the components to be processed and their beams.",
                        default='smc_srcs_image_params.vot')
    parser.add_argument("-m", "--max_loops", help="The maximum number of processing loops the daemon will run.",
                        type=int, default=500)
    parser.add_argument("-c", "--concurrency_limit", help="The maximum number of concurrent processes allowed to run.",
                        type=int, default=12)
    parser.add_argument("-n", "--min_concurrency_limit", help="The minumum number of concurrent processes we prefer to run. " +
                        "Duplicate ms usage will be allowed in order to reach this number of jobs",
                        type=int, default=6)

    parser.add_argument("--pbs", help="Run the jobs via PBS qsub command", default=False,
                        action='store_true')
    parser.add_argument("-l", "--log_folder", help="The folder which will contain the stdout and stderr files from the jobs",
                        default='logs')
    parser.add_argument("-a", "--active", help="The numerical component index of an active cutout job. The job will be monitored as if this daemon started it",
                        type=int, action='append')
    parser.add_argument("-t", "--target", help="The numerical component index of a target to be processed. If any targets are specified then only those targets are run. Default is for all targets to be run.",
                        type=int, action='append')
    parser.add_argument("--target_file", help="A file listing the names of targets to be processed. If any targets are specified then only those targets are run.",
                        type=str)
    parser.add_argument("--retry_failed", help="Cleanup any already failed jobs and retry them", default=False,
                        action='store_true')
    args = parser.parse_args()
    return args


def run_os_cmd(cmd, failOnErr=True):
    """
    Run an operating system command ensuring that it finishes successfully.
    If the comand fails, the program will exit.
    :param cmd: The command to be run
    :return: None
    """
    print(">", cmd)
    sys.stdout.flush()
    try:
        retcode = subprocess.call(cmd, shell=True)
        if retcode != 0:
            message = "Command '"+cmd+"' failed with code " + str(retcode)
            print(message, file=sys.stderr)
            if failOnErr:
                raise CommandFailedError(message)
    except OSError as e:
        message = "Command '" + cmd + "' failed " + e
        print(message, file=sys.stderr)
        if failOnErr:
            raise CommandFailedError(message)
    return None


def get_source_list(filename):
    """
    Read the sources and the beams they can be found in from the targets csv file.

    Parameters
    ----------
    filename: str
        Name of the csv file listing all sources and the beams they can be found in.
    
    Returns
    -------
    List of the sources in the same order as the csv file and a map of source names to beam lists.
    """
    src_beam_map = dict()
    targets = []
    with open(filename, 'r') as csvfile:
        tgt_reader = csv.reader(csvfile, delimiter=',',
                                quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for row in tgt_reader:
            if (tgt_reader.line_num == 1):
                # skip header
                continue
            comp_name = row[1]
            #ra = float(row[2])
            #dec = float(row[3])
            beams = row[4:]
            targets.append(comp_name)
            src_beam_map[comp_name] = beams
    return targets, src_beam_map


def read_image_params(filename):
    # targets - array with an entry per source - has a 'component_name' entry per row
    # image_params - array with component_name and beam_ids entries per row - one entry per source/beam combo
    table = parse_single_table(filename, pedantic=False)
    image_params = table.array
    targets = []
    for row in image_params:
        if row['component_name'] not in targets:
            targets.append(row['component_name'])
    return targets, image_params


def build_map(image_params):
    # Build map of sources to beam ids
    src_beam_map = dict()
    for row in image_params:
        comp_name = row['component_name']
        beam_id = row['beam_ids']
        if comp_name not in src_beam_map.keys():
            beams = set()
            src_beam_map[comp_name] = beams
        beams = src_beam_map[comp_name]
        beams.add(beam_id)
    return src_beam_map


def register_active(targets, src_beam_map, active_ids, active_ms, pre_active_jobs, remaining_array_ids, status_folder):
    if pre_active_jobs:
        active_ids.update(pre_active_jobs)

    for array_id in remaining_array_ids:
        if os.path.isfile('{}/{:d}.ACTIVE'.format(status_folder, array_id)):
            active_ids.add(array_id)

    for array_id in active_ids:
        comp_name = targets[array_id-1]
        tgt_ms = src_beam_map[comp_name]

        for ms in tgt_ms:
            active_ms.append(ms)
        # print ('+++ ' + str(active_ids))
        print('Registered active job {} (#{}) concurrency {} ms: {}'.format(
            comp_name, array_id, len(active_ids), tgt_ms))
    return len(active_ids)


def mark_comp_done(array_id, tgt_ms, active_ids, active_ms):
    active_ids.remove(array_id)
    for ms in tgt_ms:
        active_ms.remove(ms)


def job_loop(targets, sbid, status_folder, src_beam_map, active_ids, active_ms, remaining_array_ids, completed_srcs,
             failed_srcs, concurrency_limit, min_concurrency_limit, use_pbs, log_folder):
    rate_limited = False
    # Take a copy of the list to avoid issues when removing items from it
    ids_to_scan = list(remaining_array_ids)
    # Scan for completed jobs
    for array_id in ids_to_scan:
        comp_name = targets[array_id-1]
        if comp_name in completed_srcs or comp_name in failed_srcs:
            continue

        if os.path.isfile('{}/{:d}.COMPLETED'.format(status_folder, array_id)):
            # print ('--- ' + str(active_ids))
            if array_id in active_ids:
                print('Completed {}  (#{}) concurrency {}'.format(comp_name, array_id, len(active_ids)))
                mark_comp_done(array_id, src_beam_map[comp_name], active_ids, active_ms)
            else:
                print(' Skipping {} (#{}) as it has already completed'.format(comp_name, array_id))
            completed_srcs.add(comp_name)
            remaining_array_ids.remove(array_id)
            continue

        if os.path.isfile('{}/{:d}.FAILED'.format(status_folder, array_id)):
            if array_id in active_ids:
                print('Failed {}  (#{}) concurrency {}'.format(comp_name, array_id, len(active_ids)))
                mark_comp_done(array_id, src_beam_map[comp_name], active_ids, active_ms)
            else:
                print(' Skipping {} (#{}) as it has already failed'.format(comp_name, array_id))
            failed_srcs.add(comp_name)
            remaining_array_ids.remove(array_id)
            continue

    # Scan for jobs to start
    ids_to_scan = list(remaining_array_ids)
    for array_id in ids_to_scan:
        if array_id in active_ids:
            #print ('{} is active'.format(array_id))
            continue

        comp_name = targets[array_id-1]
        if comp_name in completed_srcs:
            print ('{} (#{}) has already completed as a different id!'.format(comp_name, array_id))
            remaining_array_ids.remove(array_id)
            continue


        tgt_ms = src_beam_map[comp_name]
        if len(active_ids) > min_concurrency_limit:
            clash = False
            for ms in tgt_ms:
                if ms in active_ms:
                    clash = True
            if clash:
                continue

        if len(active_ids) < concurrency_limit:
            rate_limited = False
            for ms in tgt_ms:
                active_ms.append(ms)
            active_ids.add(array_id)
            # print ('+++ ' + str(active_ids))
            print('Starting {} (#{}) concurrency {} ms: {}'.format(
                comp_name, array_id, len(active_ids), tgt_ms))
            # run_os_cmd('./make_askap_abs_cutout.sh {} {}'.format(array_id, status_folder))
            me = Path(__file__)
            script = Path(me.parent, 'start_job.sh')

            if use_pbs:
                run_os_cmd(
                    ('qsub -v COMP_INDEX={0},SBID={2},STATUS_DIR={3} -N "ASKAP_abs{0}" -o {1}/askap_abs_{0}_o.log '
                     '-e {1}/askap_abs_{0}_e.log {4}').format(array_id, log_folder, sbid, status_folder, script))
            else:
                run_os_cmd('{} {} {} "{}"'.format(script, array_id, sbid, status_folder))
        elif not rate_limited:
            rate_limited = True
            print (' rate limit of {} applied'.format(concurrency_limit))
    return len(active_ids)


def produce_all_cutouts(targets, sbid, status_folder, src_beam_map, delay, concurrency_limit, min_concurrency_limit, use_pbs,
                        log_folder, pre_active_jobs, target_list, max_loops=500):
    remaining_array_ids = list(range(1, len(targets)+1))
    if target_list:
        remaining_array_ids = list(target_list)
    num_targets = len(remaining_array_ids)
    active_ms = list()
    active_ids = set()
    completed_srcs = set()
    failed_srcs = set()

    total_concurrency = 0
    print('Processing {} targets'.format(len(remaining_array_ids)))

    num_running = register_active(targets, src_beam_map, active_ids, active_ms, pre_active_jobs, remaining_array_ids, 
                                    status_folder)
    total_concurrency += num_running
    i = 0
    while len(remaining_array_ids) > 0 and i < max_loops:
        i += 1
        print("\nLoop #{}, completed {} failed {} running {} remaining {} at {}".format(
            i, len(completed_srcs), len(failed_srcs), num_running, len(remaining_array_ids)-num_running, 
            time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))), flush=True)
        num_running = job_loop(targets, sbid, status_folder, src_beam_map, active_ids, active_ms, remaining_array_ids,
                                      completed_srcs, failed_srcs, concurrency_limit, min_concurrency_limit, use_pbs, 
                                      log_folder)
        total_concurrency += num_running
        if len(remaining_array_ids) > 0:
            time.sleep(delay)

    if len(remaining_array_ids) > 0:
        msg = 'ERROR: Failed to complete processing after {} loops. {} cutouts remain'.format(
            i, len(remaining_array_ids))
        print('\n'+msg)
        raise Exception(msg)
    else:
        print('\nCompleted processing in {} loops with average concurrency {:.2f}'.format(
            i, total_concurrency/i))
    return num_targets
    

def prep_folders(folders):
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)
            print ("Created " + folder)


def build_target_list(targets, target_ids, target_file):
    """
    Build a list of the targets to be processed.

    Parameters
    ----------
    targets: str[]
        List of the component names of possible targets.
    target_ids: int[]
        List of the numerical ids of the subset of targets to be processed.
    target_file: str
        Name of the file listing the comonent names of the subset of targets to be processed.
    
    Returns
    -------
    List of the numerical ids of the subset of targets to be processed.
    """
    if target_ids:
        return target_ids
    if not target_file:
        return None

    target_list = []
    with open(target_file, 'r') as tgt_file:
        for line in tgt_file:
            comp_name = line.strip()
            if len(comp_name) == 0 or comp_name.startswith('#'):
                continue
            found = False
            for idx, tgt in enumerate(targets):
                if tgt == comp_name:
                    target_list.append(idx+1)
                    found = True
                    continue
            if not found:
                print ("Source {} is not known. Ignoring.".format(comp_name))
    
    # Make sure we don't run everything if a list was supplied but nothing is found
    if len(target_list) == 0:
        msg = 'ERROR: None of the targets listed in {} were found.'.format(target_file)
        print (msg)
        raise Exception(msg)
    
    return target_list

def cleanup_failed(status_folder):
    failed_list = glob.glob('{}/*.FAILED'.format(status_folder))
    count = 0
    for filename in failed_list:
        id = int(os.path.basename(filename).split('.')[0])
        print ("Cleaning up failed job #{}".format(id))
        os.remove(filename)
        count+= 1
    print ("Cleaned up {} failed jobs.".format(count))


def main():
    # Parse command line options
    args = parseargs()

    start = time.time()
    print("#### Started ASKAP cutout production of sbid {} at {} ####".format
          (args.sbid, time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start))))

    print ('Checking every {} seconds for completed jobs, with a maximum of {} checks.'.format(args.delay, args.max_loops))

    work_folder = 'sb{}/work'.format(args.sbid)
    cutouts_folder = 'sb{}/cutouts'.format(args.sbid)

    print (' Status folder', args.status_folder)
    print (' Log folder', args.log_folder)
    print (' Source filename', args.filename)
    print (' Target list', args.target)
    print (' Concurrency max {} min {}'.format(args.concurrency_limit, args.min_concurrency_limit))
    print (' Batch system', ('PBS' if args.pbs else 'None'))

    # Prepare the run
    #targets, image_params = read_image_params(args.filename)
    #src_beam_map = build_map(image_params)
    targets, src_beam_map = get_source_list(args.filename)
    target_list = build_target_list(targets, args.target, args.target_file)

    if args.active:
        print ("\nAlready active jobs: {}".format(args.active))
    if target_list:
        print ("\nLimiting run to only these targets/jobs: {}".format(target_list))
    
    status_folder = '{}/{}'.format(args.status_folder, args.sbid)
    log_folder = '{}/{}'.format(args.log_folder, args.sbid)
    prep_folders([status_folder, log_folder, work_folder, cutouts_folder])

    if args.retry_failed:
        cleanup_failed(status_folder)

    # Run through the processing
    num_targets = produce_all_cutouts(targets, args.sbid, status_folder, src_beam_map, args.delay, 
                        args.concurrency_limit, args.min_concurrency_limit, args.pbs, 
                        log_folder, args.active, target_list, max_loops=args.max_loops)

    # Report
    end = time.time()
    print('#### Processing completed at {} ####'.format(
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end))))
    print('Processed {0} targets in {1:.2f} min'.format(
          num_targets, (end - start)/60))
    return 0


if __name__ == '__main__':
    exit(main())
