#!/usr/bin/env python -u

# Daemon process to manage the production of cutouts from an ASKAP scheduing block

# Author James Dempsey
# Date 27 Jan 2020

import argparse
import datetime
# import math
import os
import subprocess
import sys
import time
# from collections import deque

# from astropy.io import ascii
# from astropy.coordinates import SkyCoord, Angle
# from astropy.io import fits, votable
# from astropy.wcs import WCS
# import numpy as np
# import numpy.core.records as rec
# from astropy import units as u
# from astropy.table import QTable, Table, Column
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
    parser.add_argument("-s", "--status_folder", help="The status folder which will contain the completed files",
                        default='status/8906')
    parser.add_argument("-f", "--filename", help="The name of the votabel format file listing the components to be processed.",
                        default='smc_srcs_image_params.vot')
    parser.add_argument("-m", "--max_loops", help="The maximum number of processing loops the daemon will run.",
                        type=int, default=500)
    parser.add_argument("-c", "--concurrency_limit", help="The maximum number of concurrent processes allowed to run.",
                        type=int, default=12)
    parser.add_argument("--pbs", help="Run the jobs via PBS qsub command", default=False,
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


def job_loop(targets, status_folder, src_beam_map, active_ids, active_ms, remaining_array_ids, completed_srcs, concurrency_limit, use_pbs):
    rate_limited = False
    # Take a copy of the list to avoid issues when removing items from it
    ids_to_scan = list(remaining_array_ids)
    # Scan for completed jobs
    for array_id in ids_to_scan:
        comp_name = targets[array_id-1]
        if comp_name in completed_srcs:
            continue

        if os.path.isfile('{}/{:d}.COMPLETED'.format(status_folder, array_id)):
            # print ('--- ' + str(active_ids))
            if array_id in active_ids:
                print('Completed {}  (#{}) concurrency {}'.format(
                    comp_name, array_id, len(active_ids)))
                active_ids.remove(array_id)
                tgt_ms = src_beam_map[comp_name]
                for ms in tgt_ms:
                    active_ms.remove(ms)
            else:
                print(' Skipping {} (#{}) as it has already completed'.format(
                    comp_name, array_id))
            completed_srcs.add(comp_name)
            remaining_array_ids.remove(array_id)
            continue

    # Scan for jobs to start
    for array_id in remaining_array_ids:
        comp_name = targets[array_id-1]
        if comp_name in completed_srcs:
            continue

        tgt_ms = src_beam_map[comp_name]
        if tgt_ms & active_ms:
            continue

        if len(active_ids) < concurrency_limit:
            rate_limited = False
            for ms in tgt_ms:
                active_ms.add(ms)
            active_ids.add(array_id)
            # print ('+++ ' + str(active_ids))
            print('Starting {} (#{}) concurrency {} ms: {}'.format(
                comp_name, array_id, len(active_ids), tgt_ms))
            # run_os_cmd('./make_askap_abs_cutout.sh {} {}'.format(array_id, status_folder))
            if use_pbs:
                run_os_cmd('qsub -J {0}-{1}:2 -N "ASKAP_abs{0}" ./start_job.sh'.format(array_id,array_id+1))
            else:
                run_os_cmd('./start_job.sh {}'.format(array_id))
        elif not rate_limited:
            rate_limited = True
            print (' rate limit of {} applied'.format(concurrency_limit))
    return len(active_ids)


def produce_all_cutouts(targets, status_folder, src_beam_map, delay, concurrency_limit, use_pbs, max_loops=500):
    remaining_array_ids = list(range(1, len(targets)+1))
    active_ms = set()
    active_ids = set()
    completed_srcs = set()

    total_concurrency = 0
    print('Processing {} targets'.format(len(remaining_array_ids)))
    i = 0
    while len(remaining_array_ids) > 0 and i < max_loops:
        i += 1
        print("\nLoop #{}, completed {} remaining {} at {}".format(
            i, len(completed_srcs), len(remaining_array_ids), time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(time.time()))), flush=True)
        total_concurrency += job_loop(targets, status_folder, src_beam_map, active_ids, active_ms, remaining_array_ids,
                                      completed_srcs, concurrency_limit, use_pbs)
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


def main():
    # Parse command line options
    args = parseargs()

    start = time.time()
    print("#### Started ASKAP cutout production at %s ####" %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(start)))

    # Prepare the run
    targets, image_params = read_image_params(args.filename)
    src_beam_map = build_map(image_params)

    # Run through the processing
    produce_all_cutouts(targets, args.status_folder,
                        src_beam_map, args.delay, args.concurrency_limit, args.pbs, max_loops=args.max_loops)

    # Report
    end = time.time()
    print('#### Processing completed at %s ####' %
          time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(end)))
    print('Processed %d components in %.02f s' %
          (len(targets), end - start))
    return 0


if __name__ == '__main__':
    exit(main())
