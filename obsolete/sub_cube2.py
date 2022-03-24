vis = []
for beam_offset in (['07', 'A'], ['02', 'B'], ['07', 'C']):
    vis_name='scienceData_SB8906_SMC1-0_M344-11{}.beam{}_SL.ms'.format(beam_offset[1], beam_offset[0])
    vis.append(vis_name)

print "Processing inpout visibilities of " + str(vis)

for klambda in ('0.8', '0.6', '0.4', '0.2'):
        print ('uvdist={}Klambda'.format(klambda))
	image_name='sb8906/sb8906_0029-7228_sl_{}K'.format(klambda)
	uvdist='>{}Klambda'.format(klambda)
	fits_name='sb8906/sb8906_0029-7228_sl_{}K.fits'.format(klambda)
	tclean (vis=vis,specmode='cube',imagename=image_name,reffreq='1.42040571183GHz',restfreq='1.42040571183GHz',
	  phasecenter='J2000 00h29m19 -72d28m11',imsize=50,uvrange=uvdist,
	  gridder='standard', width='1km/s',
	  vptable='ASKAP_AIRY_BP.tab')
	exportfits(imagename=image_name+'.image', fitsimage=fits_name,velocity=True)

