give_me_info_lim.pro
;    -Convert UTC to JULIAN to MY
;    -Get attributes from respective shapefiles
;    -Comments in all caps describe actions
;


PRO give_me_info


;==== SELECT DETECTOR TYPE AND OBS MODE ====;
dtype = ''
 read, dtype, prompt = 'Enter detector type: '
  
obs_mode = ''
 read, obs_mode, prompt = 'Enter observing mode (ALL CAPS): '


;==== Target file root ====;
target_dbf_file_root_root = '/disks/dc026/msp_mos/north_polar_mapping/'
target_dbf_file_root = target_dbf_file_root_root + dtype
  print, 'BEGIN SEARCH FOR FILES: ' + target_dbf_file_root + '/EVF/' + obs_mode + '/*.dbf'

target_dbf_file = file_search(target_dbf_file_root + '/EVF/' + obs_mode + '/*.dbf', count = count1)
  print, count1


;==== PULL ATTRIBUTES FROM THE SHAPEFILE ====;
crism_dbf = obj_new('IDLffShape', target_dbf_file, /dbf_only)
crism_dbf->GetProperty, ATTRIBUTE_NAMES = attribute_names, ATTRIBUTE_INFO = attribute_info
attributes = crism_dbf->IDLffShape::getattributes(/all)
 ;print, attributes

  ;stop
  
num_elements_in_file = n_elements(attributes)
 print, num_elements_in_file


;==== Specific attribute reference list ====;
A1 = attributes.(0); A0_filename
A2 = attributes.(5); A5_classtype
A3 = attributes.(6); A6_obsID -- targetID
A4 = attributes.(7); A7_counter -- segmentID
A5 = attributes.(19); A19_starttime
A6 = attributes.(24); A24_solarlongi
A7 = attributes.(27); A27_detector

filename = attributes.(0); A0_filename
classtype = attributes.(5); A5_classtype
targetID = attributes.(6); A6_obsID -- targetID
segmentID = attributes.(7); A7_counter -- segmentID
starttime_utc = attributes.(19); A19_starttime
slongi = attributes.(24); A24_solarlongi
detector = attributes.(27); A27_detector


;==== WRITE OUT ATTRIBUTES and here's a reference list too ====;

    write_my_attributes, dtype = dtype, obs_mode = obs_mode, A1 = A1, A2 = A2, A3 = A3, A4 = A4, A5 = A5, A6 = A6, A7 = A7, $
                         filename = filename, classtype = classtype, targetID = targetID, segmentID = segmentID, $
                         starttime = starttime, slongi = slongi, detector = detector

working_directory = '/disks/dc026/msp_mos/north_polar_mapping/' + dtype + '/SLF/PRJ/' + obs_mode + '_attributes/'
  print, 'Working with attribute files located in: ' + working_directory
 
A1_list_file = working_directory + 'filenames.txt'
A2_list_file = working_directory + 'classtypes.txt'
A3_list_file = working_directory + 'targetIDs.txt'
A4_list_file = working_directory + 'segmentIDs.txt'
A5_list_file = working_directory + 'starttimes.txt'
A6_list_file = working_directory + 'slongis.txt'
A7_list_file = working_directory + 'detector.txt'

; *_FIXED files are corrected for the DOY default reading in mro_crism_utc_julian_convert (happened in the 1st run of give_me_info_really)
A5_Julian_list_file_fixed = working_directory + 'starttimes_julian_fixed.txt' ;A5 in Julian
A5_MY_list_file_fixed = working_directory + 'starttimes_my_fixed.txt' ;A5 in MY

; *_REALLY files are corrected for the row of data (need a column instead)
A5_Julian_list_file_really = working_directory + 'starttimes_julian_really.txt' ;A5 in Julian
A5_MY_list_file_really = working_directory + 'starttimes_MY_really.txt' ;A5 in MY

; IF FILES DON'T EXIST: ******* write_my_attributes
;THEN you need to do this to format it into 1 column of info
;cat A*_list.txt | xargs -n 1 > A*_list.txt
;sed '/^$/d' A*_list.txt > A*_list.txt


;stop

  obj_destroy, crism_dbf


;=============== UTC TO JULIAN ================

Ju_file_search = file_search(A5_Julian_list_file_fixed, count = num_Ju_files_fixed)

if num_Ju_files_fixed EQ 1 then begin
  print, 'Start times have already been converted from UTC to Julian (in ymd format)'

endif else begin
            
  print, 'GIVE_ME_INFO_REALLY: Begin reading in UTC files...'
            
  openr, UTC_lun, A5_list_file, /get_lun
  readf, UTC_lun, UTC
  free_lun, UTC_lun
             
  print, 'Done :)'           
  print, 'GIVE_ME_INFO_REALLY: Begin converting from UTC to Julian...'
               
UTC_to_Ju = mro_crism_utc_julian_convert(utc = A5, /ymd)

  openw, Jul_lun, A5_Julian_list_file_fixed, /get_lun
  printf, Jul_lun, n_elements(A5)
  printf, Jul_lun, UTC_to_Ju.julian, format = '(f12.4)'
  free_lun, Jul_lun
               
  print, 'Done :)'
  
      endelse  

 
;================ JULIAN TO MY ==================

MY_file_search = file_search(A5_MY_list_file_really, count = num_MY_files)
 
 if num_MY_files EQ 1 then begin
  print, 'List A5 has already been converted from Julian to MY'
 endif else begin
 
  if num_MY_files LT 1 then begin         
   print, 'Begin reading in Julian files...'          
    readcol, A5_Julian_list_file_fixed, Julian, skipline = 1
    
;stop

     print, 'Begin converting from Julian to MY...'        
      num_Julian_dates = n_elements(Julian)
      Ju_to_MY = make_array(num_Julian_dates, value = 0)
        
        for i = 0, num_Julian_dates -1 do begin
         Ju_to_MY[i] = mro_crism_mars_year(julian_date = Julian[i])
        endfor
     
     openw, MY_lun, A5_MY_list_file_really, /get_lun
     printf, MY_lun, n_elements(A5)
     printf, MY_lun, transpose(Ju_to_MY)
     free_lun, MY_lun
     print, 'Done :)'
 
            
  endif
 
 
        
 endelse 

    parse_my_info, dtype = dtype, obs_mode = obs_mode

END                                                                                                                                                                                                                                                                                                               mosaic_root_and_name.txt                                                                            0000664 €    80002301 00000000151 13120770765 015257  0                                                                                                    ustar   limh1                           crism-soc                                                                                                                                                                                                              mosaic_root = '/project/crism/users/seelofp1/msp_mos/north_polar_mapping/'

mosiac_name = 'planum_boreum'
                                                                                                                                                                                                                                                                                                                                                                                                                       mro_crism_build_virtual_mosaic_lim.pro                                                              0000664 €    80002301 00000014054 13144367416 020223  0                                                                                                    ustar   limh1                           crism-soc                                                                                                                                                                                                              ;06/27/17 (limh1)
;   -Added keywords for detector and processing paths
;   -Added program to write the mosaic files in envi format
;   -This is particularly for the parsed out 10degree L_s windows, either full or partial

pro mro_crism_build_virtual_mosaic_lim, part_mos = part_mos, quick ;, full_mos = full_mos

;To compile envi functions w/o having envi running
compile_opt strictarr


dtype = ''
  read, dtype, prompt = 'Enter detector type (ir or vnir): '
    
    if dtype EQ 'ir' then begin
      dtype_short = 'L'
        endif else begin
          dtype_short = 'S'
        endelse

obs_mode = ''
  read, obs_mode, prompt = 'Enter observing mode (HSP or MSP): '


if keyword_set(quick) then begin
  C_list_file = ''
   read, C_list_file, prompt = 'Enter full path to list file: '
endif


if keyword_set(part_mos) then begin 
 Lower_longi = ''
  read, Lower_longi, prompt = 'Enter lower longitude bound (ie. 150): '
 Upper_longi = ''
  read, Upper_longi, prompt = 'Enter upper longitude bound (ie. 160): '
 
  print, 'MRO_CRISM_BUILD_VIRTUAL_MOSAIC_LIM: Begin building for Solar Longitude' + Lower_longi + ' to ' + Upper_longi + '...'
 
 C_list_file = '/disks/dc026/msp_mos/north_polar_mapping/' + dtype + '/SLF/PRJ/' + obs_mode + '_attributes/Ls/Ls_' + Lower_longi + '_' + Upper_longi + '.txt'
 
endif


 C_list_struct = mro_crism_read_list(C_list_file)
  
  stop
 
 DDR_list_file = '/disks/dc026/msp_mos/north_polar_mapping/' + dtype + '/EVF/' + obs_mode + '/planum_boreum_' + obs_mode + '_' + dtype_short + '_EVF_INT.txt'
 DDR_list_struct = mro_crism_read_list(DDR_list_file)
 DDR_file_struct = mro_crism_filename_struct_array(DDR_list_struct.(1))
 DDR_file_struct_ident = mro_crism_parse_identifier(filename_struct = DDR_file_struct)
 DDR_indx = cmset_op(DDR_file_struct_ident, 'AND', C_list_struct.(1), /index)
  
 DDR_input_file_list = DDR_list_struct.(1)[DDR_indx]
 ddr_file_list = DDR_input_file_list
  
 PHT_list_file = '/disks/dc026/msp_mos/north_polar_mapping/' + dtype + '/PHT/' + obs_mode + '/planum_boreum_' + obs_mode + '_' + dtype_short + '_TRR_PRE_PHT.txt'
 PHT_list_struct = mro_crism_read_list(PHT_list_file)
 PHT_file_struct = mro_crism_filename_struct_array(PHT_list_struct.(1))
 PHT_file_struct_ident = mro_crism_parse_identifier(filename_struct = PHT_file_struct)
 PHT_indx = cmset_op(PHT_file_struct_ident, 'AND', C_list_struct.(1), /index)
 PHT_input_file_list = PHT_list_struct.(1)[PHT_indx]
 trdr_file_list = pht_input_file_list
  

 ddr_output_files = '/disks/dc026/msp_mos/north_polar_mapping/' + dtype + '/SLF/PRJ/xovr_' + Lower_longi + '_' + Upper_longi + '_' + obs_mode + '_DDR.txt'
  openw, dlun, ddr_output_files, /get_lun
  printf, dlun, DDR_input_file_list
  free_lun, dlun

 pht_output_file = '/disks/dc026/msp_mos/north_polar_mapping/' + dtype + '/SLF/PRJ/xovr_' + Lower_longi + '_' + Upper_longi + '_' + obs_mode + '_PRE_PHT.txt'
  openw, plun, pht_output_file, /get_lun
  printf, plun, PHT_input_file_list
  free_lun, plun
  
 input_list_file = pht_output_file
  ;print, 'MRO_CRISM_BUILD_VIRTUAL_MOSAIC_LIM: Building virtual mosaics for ' + input_files_count + 'files'


envi_check = envi(/current)
if (obj_valid(envi_check) EQ 0) then begin
  print, 'Initalizating headless ENVI session...'
    envi_object = envi(/headless)
  print, 'Done!'
    ;stop
endif


;input_file_list = file_search('/disks/dc026/msp_mos/north_polar_mapping/' + detector_path + 'SLF/PRJ/' + processing_path + '*.IMG', count=count1)
;input_list_file = file_search('/disks/dc026/msp_mos/north_polar_mapping/' + dtype + '/SLF/PRJ/xovr_*.txt', count = input_files_count)
;print, 'MRO_CRISM_BUILD_VIRTUAL_MOSAIC_LIM: Building virtual mosaics for ' + input_files_count + 'files'

    readcol, input_list_file, xover_files, format = 'A'
 num_input_list = n_elements(xover_files)
  
stop

 ;num_input_list = n_elements(input_list_file)]
 input_fid_list = make_array(num_input_list, value = -1, /long)
 input_dims_array = make_array(5, num_input_list, value = -1, /long)
  
 for i = 0, num_input_list -1 do begin
  
   ;envi_open_file, input_file_list[i], r_fid = input_fid, /no_interactive_query, /no_realize
   envi_open_file, xover_files[i], r_fid = input_fid, /no_interactive_query, /no_realize
   envi_file_query, input_fid, dims = input_dims, ns = input_ns, nl = input_nl, nb = input_nb, fname = input_filename
    
  input_dims_array[*,i] = input_dims
  input_fid_list[i] = input_fid
  ;Retrieve band names and projection info from first input file
   envi_file_query, input_fid, bnames = input_bnames
   print, input_bnames

  input_proj = envi_get_projection(fid = input_fid, pixel_size = input_pxl_size)
   help, input_proj, /struct
	 print, input_pxl_size
 
stop
       
 endfor


georef_mosaic_setup, fids = input_fid_list, dims = input_dims_array, out_ps = input_pxl_size, $
    	    	         xsize = mosaic_xsize, ysize = mosaic_ysize, x0 = mosaic_x0, y0 = mosaic_y0, $
                     map_info = mosaic_map_info
                    
see_through_val = make_array(n_elements(xover_files), value = 65535.)
;print, see_through_val

;stop

write_mosaic_file_lim, output_file = output_file, input_file_list = input_file_list, input_list_file = input_list_file, $
                      input_fid_list = input_fid_list, input_dims = input_dims_array, $
                      input_bnames = input_bnames, mosaic_xsize = mosaic_xsize, $
                      mosaic_ysize = mosaic_ysize, input_pxl_size = input_pxl_size, $
                      mosaic_x0 = mosaic_x0, mosaic_y0 = mosaic_y0, see_through_val = see_through_val, xover_files = xover_files
                   

;envi_doit, 'mosaic_doit', fid = input_fid_list, pos = pos, dims = input_dims_array, out_name = out_name, $
 ;                  xsize = mosaic_xsize, ysize = mosaic_ysize, x0 = mosaic_x0, y0 = mosaic_y0, georef = 1, map_info = mosaic_map_info, $
  ;                 out_dt = data_type, pixel_size = input_pxl_size, background = 65535., see_through_val = VSTV, $
   ;                ;use_see_through = use_see_through,
    ;               out_bname = input_bnames
 
if (obj_valid(envi_check) EQ 0) then begin
    print, 'Closing headless ENVI session...'
    envi_object.close
endif



END
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    mro_crism_glt_generate_lim.pro                                                                      0000664 0003330 0002301 00000025474 13122533300 017345  0                                                                                                    ustar   seelofp1                        crism-soc                                                                                                                                                                                                              ;10/17/2007 (fps)
;   Updated for use in revamped sandbox mosaicking environment
;   Keyword geo_subset changed from logical to lat/lon bounds array (if active)
;   Parsing of geo_subset array ensures scalar constituent variables
;   Disabled list guided product update (for now...)
;   Removed output_path and output_list_file redundancy - output_list_file is written to the output path
;   Improved trap for single spatial pixel output domain
;

pro mro_crism_glt_generate_lim, output_path, in_proj, out_proj, input_list_file = input_list_file, input_file_list = input_file_list, output_list_file = output_list_file, product_check = product_check, $
    	    	    	    compression = compression, geo_subset = geo_subset, trim_subset = trim_subset, pxl_size = pxl_size
			    ;lon_min = lon_min, lon_max = lon_max, lat_min = lat_min, lat_max = lat_max

compile_opt strictarr

help, in_proj, /struct
help, out_proj, /struct

if (keyword_set(input_file_list)) then begin
    num_input_ddr_list = n_elements(input_file_list)
    input_ddr_list = input_file_list
endif


if (keyword_set(input_list_file)) then begin

    input_list_file_info = file_info(input_list_file)

    if (input_list_file_info.exists EQ 1) then begin
    	input_list_file_struct = mro_crism_msp_mosaic_read_list(input_list_file)
    	num_input_ddr_list = (input_list_file_struct).(0)
    	input_ddr_list = (input_list_file_struct).(1)
    endif
    
endif

;Trying to split each HSP, MSP and MSW .txt file into .txt files with size > or = 500 lines

;if (keyword_set(split_list)) then begin
;    num_input_ddr_list = n_elements(split_list)
;    input_ddr_list = split_list
;endif


input_ddr_file_struct_array = mro_crism_filename_struct_array(input_ddr_list)


num_extant_glt_list = 0
;extant_glt_list = ''
;
;if (keyword_set(output_list_file)) then begin
;
;    output_list_file_info = file_info(output_list_file)
;
;    if (output_list_file_info.exists EQ 1) then begin
;    	output_list_file_struct = mro_crism_msp_mosaic_read_list(output_list_file)
;    	num_extant_glt_list = (output_list_file_struct).(0)
;    	extant_glt_list = (output_list_file_struct).(1)
;    endif
;
;endif

print, output_path
print, input_ddr_list
;print, extant_glt_list

help, input_ddr_file_struct_array, /struct


if (num_extant_glt_list GT 0) then begin
    extant_glt_file_struct_array = mro_crism_filename_struct_array(extant_glt_list)
    
    for i = 0, n_elements(extant_glt_file_struct_array) -1 do begin
    
    endfor
    
endif else begin

    ddr_glt_array = input_ddr_list
    num_ddr_glt_array = num_input_ddr_list

    active_glt_array = make_array(num_ddr_glt_array, value = '')
    active_glt_array_indx = replicate(long(1), num_ddr_glt_array)

endelse

envi_check = envi(/current)
		if (obj_valid(envi_check) EQ 0) then begin
			print, 'MRO_CRISM_GLT_GENERATE_UPDATE: Initalizating headless ENVI session...'
		envi_object = envi(/headless)
			print, 'Done!'
  	;stop
		endif


for i = 0, num_ddr_glt_array -1 do begin

    print, 'Generating GLT ' + strtrim(string(i+1),2) + ' / ' + strtrim(string(fix(num_ddr_glt_array)),2) + '     ' + systime(0)

    ddr_filename = ddr_glt_array[i]
    ddr_filename_struct = mro_crism_parse_filename(ddr_filename)

    glt_filename =  output_path + ddr_filename_struct.class_type + ddr_filename_struct.obs_id + '_' + ddr_filename_struct.counter + '_' + $
    	    	    ddr_filename_struct.activity + ddr_filename_struct.sensor + '_' + 'GLT' + ddr_filename_struct.version + ddr_filename_struct.file_extension

    print, glt_filename
    active_glt_array[i] = glt_filename

    glt_run_flag = 1
    if (keyword_set(product_check)) then begin
    	glt_run_flag = mro_crism_overwrite_check(glt_filename)
    endif 

    ;stop

    if (glt_run_flag EQ 1) then begin

	ddr = mro_crism_quick_read(ddr_filename)
	ddr_lon = ddr[*,*,4]
	ddr_lat = ddr[*,*,3]

	ddr_size = size(ddr)
	print, ddr_size

	;Enter DDR data (lat, lon, etc.) into ENVI ABL
	envi_enter_data, ddr, r_fid = ddr_fid

	;Query newly entered file for ENVI info
	envi_file_query, ddr_fid, dims = ddr_dims, ns = ddr_ns, nl = ddr_nl, nb = ddr_nb

	print, ddr_dims
	print, ddr_ns, ddr_nl, ddr_nb

	;if (0 EQ 1) then begin
	;if (keyword_set(geo_subset)) then begin
	if (n_elements(geo_subset) EQ 4) then begin
    	    ;assume geo bounds are valid here - check further up at some point...
	    ;recover geographic bounds...
	    lon_min = geo_subset[0]
	    lon_max = geo_subset[1]
	    lat_min = geo_subset[2]
	    lat_max = geo_subset[3]

       	    ;Modify ddr_dims to subset DDR specified lat/lon region for GLT construction

	        ddr_lon_active_indx = where((ddr_lon GE lon_min) AND (ddr_lon LE lon_max), num_ddr_lon_active_indx)	
    	    ddr_lat_active_indx = where((ddr_lat GE lat_min) AND (ddr_lat LE lat_max), num_ddr_lat_active_indx)

    	    ddr_lon_active_mask = make_array(n_elements(ddr_lon[*,0]), n_elements(ddr_lon[0,*]), value = long(0))
    	    if (num_ddr_lon_active_indx GT 0) then begin
    		print, 'NUM DDR LON ACTIVE INDX: ', num_ddr_lon_active_indx, ' / ', n_elements(ddr_lon_active_mask)
        	ddr_lon_active_mask[ddr_lon_active_indx] = long(1)
    	    endif

	        ddr_lat_active_mask = make_array(n_elements(ddr_lat[*,0]), n_elements(ddr_lat[0,*]), value = long(0))
    	    if (num_ddr_lat_active_indx GT 0) then begin
    		print, 'NUM DDR LAT ACTIVE INDX: ', num_ddr_lat_active_indx, ' / ', n_elements(ddr_lat_active_mask)
        	ddr_lat_active_mask[ddr_lat_active_indx] = long(1)
    	    endif

          ddr_indx = where((ddr_lat_active_mask EQ 1) AND (ddr_lon_active_mask EQ 1), num_ddr_indx)

    	    ddr_indx_array = array_indices(ddr_lon_active_mask, ddr_indx)

    	    ddr_dims = [-1, (0 > (min(ddr_indx_array[0,*]) - 1)), ((ddr_size[1] - 1) < (max(ddr_indx_array[0,*]) + 1)), $
    	    	    	    (0 > (min(ddr_indx_array[1,*]) - 1)), ((ddr_size[2] - 1) < (max(ddr_indx_array[1,*]) + 1))]

    	    ;stop

	endif
	;endif

	glt_rotation = 0.0
	;glt_pixel_size = out_proj.params[4]
	glt_pixel_size = pxl_size

    ;    stop
	if ((keyword_set(compression)) OR (keyword_set(trim_subset))) then begin
    	    ;stop
	    
            envi_doit, 'envi_glt_doit', dims = ddr_dims, i_proj = in_proj, o_proj = out_proj, $
    	    	    	    		x_fid = ddr_fid, x_pos = 4, y_fid = ddr_fid, y_pos = 3, $
					rotation = glt_rotation, pixel_size = glt_pixel_size, r_fid = glt_fid, /in_memory

    	    if (glt_fid NE -1) then begin

    		envi_file_query, glt_fid, bnames = glt_bnames, wl = glt_wl, nl = glt_nl, ns = glt_ns, nb = glt_nb, dims = glt_dims, data_type = glt_data_type, $
    	    	    		 offset = glt_offset, data_ignore_value = glt_data_ignore_value

    		glt_map = envi_get_map_info(fid = glt_fid)

    		print, 'GLT FID: ', glt_fid
		print, glt_ns, glt_nl, glt_nb

    	    	;Watch out for single line/sample ranges... very bad things..
    	    	if (total([glt_ns, glt_nl, glt_nb] EQ 1) EQ 3) then begin
		    print, 'GLT FUBAR!
    	    	    envi_file_mng, id = glt_fid, /remove
		    print, 'SETTING GLT FID = -1'
		    glt_fid = -1
		endif else begin	

    		    glt_cube = make_array(glt_ns, glt_nl, glt_nb, value = 0, type = glt_data_type)

		    for k = 0, glt_nb -1 do begin
	    ;	    print, k

    			glt_band = envi_get_data(fid = glt_fid, dims = glt_dims, pos = k)
			glt_cube[*,*,k] = glt_band
		    endfor


    		    if (keyword_set(trim_subset)) then begin

    			;generate geographic projection structure for use later
    			geo_proj_struct = envi_proj_create(/geographic)

    			;snag GLT projection info
    			glt_proj_struct = envi_get_projection(fid = glt_fid, pixel_size = glt_pxl_size, units = glt_units)    

			;calculate MAP coordinates
			glt_pxl_coords = array_indices(glt_band, lindgen(n_elements(glt_band)))
			envi_convert_file_coordinates, glt_fid, glt_pxl_coords[0,*], glt_pxl_coords[1,*], x_map, y_map, /to_map

			;calculate GEO coordinates
			envi_convert_projection_coordinates, x_map, y_map, glt_proj_struct, x_geo, y_geo, geo_proj_struct

	;    	    	stop

    			glt_trim_indx = where((x_geo LT lon_min) OR (x_geo GT lon_max) OR (y_geo LT lat_min) OR (y_geo GT lat_max), num_glt_trim_indx)

			if (num_glt_trim_indx GT 0) then begin

	    		    for k = 0, glt_nb -1 do begin 
    	    			glt_band = glt_cube[*,*,k]
    	    			glt_band[glt_trim_indx] = 0
	    			glt_cube[*,*,k] = glt_band
    	    		    endfor
    			    ;stop

			endif

		    endif


    		    envi_write_envi_file, glt_cube, ns = glt_ns, nl = glt_nl, nb = glt_nb, data_type = glt_data_type, offset = glt_offset, $
	    	    	    		      bnames = glt_bnames, wl = glt_wl, map_info = glt_map, data_ignore_value = glt_data_ignore_value, $
					      r_fid = glt_comp_fid, out_name = glt_filename, /compression, /no_open
    	    	endelse

    	    endif else begin
		print, 'GLT FID = -1'
		;stop
	    endelse

	endif else begin

            envi_doit, 'envi_glt_doit', dims = ddr_dims, i_proj = in_proj, o_proj = out_proj, out_name = glt_filename, $
    	    	    	    		x_fid = ddr_fid, x_pos = 4, y_fid = ddr_fid, y_pos = 3, $
					rotation = glt_rotation, pixel_size = glt_pixel_size, r_fid = glt_fid

	endelse

	if (glt_fid EQ -1) then begin
    	    active_glt_array[i] = ''
	    active_glt_array_indx[i] = 0
	    print, 'REJECTED...'
	    ;stop
	endif else begin
    	    envi_file_mng, id = glt_fid, /remove
	endelse

	if (ddr_fid NE -1) then begin
            envi_file_mng, id = ddr_fid, /remove
	endif

    endif   ;glt_run_flag
    
    	if ((keyword_set(select_indx))) then begin
    	    select_indx = lindgen(num_input_list)
    	endif
    	    num_select_indx = n_elements(select_indx)
        		
endfor

 if (obj_valid(envi_check) EQ 0) then begin
        print, 'MRO_CRISM_GLT_GENERATE_LIM_UPDATE: Closing headless ENVI session...'
        envi_object.close
 endif

active_glt_array = active_glt_array[where(active_glt_array_indx EQ 1)]
glt_write_flag = mro_crism_write_list(active_glt_array, output_path + output_list_file)


END

;make damn sure the lat/lon bounds are scalars
;lon_min_size = size(lon_min, /struct)
;lon_max_size = size(lon_max, /struct)
;lat_min_size = size(lat_min, /struct)
;lat_max_size = size(lat_max, /struct)
;
;if (total([lon_min_size.n_elements, lon_max_size.n_elements, lat_min_size.n_elements, lat_max_size.n_elements]) NE 4) then begin
;    print, 'Uh....no....
;    stop
;endif 
;
;if (total([lon_min_size.n_dimensions, lon_max_size.n_dimensions, lat_min_size.n_dimensions, lat_max_size.n_dimensions]) NE 0) then begin
;    print, 'Reducing dimension of input bounds...
;    lon_min = lon_min[0]
;    lon_max = lon_max[0]
;    lat_min = lat_min[0]
;    lat_max = lat_max[0]
;endif
                                                                                                                                                                                                    mro_crism_proximal_optimization_init_lim.pro                                                        0000664 €    80002301 00000042357 13144370412 021504  0                                                                                                    ustar   limh1                           crism-soc                                                                                                                                                                                                              ;07/31/2017 (limh1)
;   -Added Ls window options for partial mosaic for an appropriately named output file
;



pro mro_crism_proximal_optimization_init_lim, ddr_file_list = ddr_file_list, trdr_file_list = trdr_file_list, ddr_list_file = ddr_list_file, trdr_list_file = trdr_list_file, $
    	    	    	    	    	  restore_initial = restore_initial, restore_vectors = restore_vectors, restore_matricies = restore_matricies, part_mos = part_mos;, full_mos = full_mos

;var_file = file_search('*/var_list_for_prox_opt_init.txt')
;openr, var_lun, var_file, /get_lun
;readf, var_lun

;stop

dtype = ''
read, dtype, prompt = 'Enter detector type (ir or vnir): '

if dtype EQ 'ir' then begin
   dtype_short = 'L'
endif
if dtype EQ 'vnir' then begin
   dtype_short = 'S'
endif

obs_mode = ''
read, obs_mode, prompt = 'Enter observation mode (MSP or HSP): '

;Band and index info for matrix stuff
band = ''
read, band, prompt = 'Enter wavelength in nm (i.e. 770): '



root = '/disks/dc026/msp_mos/north_polar_mapping/'
;save_path = '/project/crism/users/seelofp1/working/'
;save_path = '/disks/dc035/sandbox/jpl_rfp/T0870/ir/BAL/'
;save_path = '/disks/dc026/msp_mos/north_polar_mapping/vnir/LIM/MSP'
save_path = root + dtype + '/LIM/' + obs_mode + '/'


if keyword_set(part_mos) then begin

  print, 'MRO_CRISM_PROXIMAL_OPTIMIZATION_INIT_LIM: Begin working on one 10 degree L_s window...'

  lower_longi = ''
  read, lower_longi, prompt = 'Enter lower longitude bound (i.e. 180): '

  ;upper_longi = lower_longi + 10
  upper_longi = ''
  read, upper_longi, prompt = 'Enter upper longitude bound (i.e. lower longitude + 10): '

  C_list_file = root + dtype + '/SLF/PRJ/' + obs_mode + '_attributes/C_Ls_' + lower_longi + '_' + upper_longi + '.txt'
  
  initial_save = save_path + 'proximal_initial_' + lower_longi + '_' + upper_longi + '.sav'
  
endif


;if keyword_set(full_mos) then begin
;  print, 'MRO_CRISM_PROXIMAL_OPTIMIZATION_INIT_LIM: Begin working on all L_s windows...' 
;  C_list_file = file_search('/disks/dc026/msp_mos/north_polar_mapping/' + dtype + '/SLF/PRJ/' + obs_mode + '_attributes/C_Ls_*.txt', count = num_C_list_files)
;initial_save = save_path + 'proximal_initial_current.sav'



C_list_struct = mro_crism_read_list(C_list_file)

;DDR_list_file = '/disks/dc026/msp_mos/north_polar_mapping/' + dtype + '/EVF/' + obs_mode + '/planum_boreum_' + obs_mode + '_' + dtype_short + '_EVF_INT.txt'
DDR_list_file = root + dtype + '/EVF/' + obs_mode + '/planum_boreum_' + obs_mode + '_' + dtype_short + '_EVF_INT.txt'
DDR_list_struct = mro_crism_read_list(DDR_list_file)
DDR_file_struct = mro_crism_filename_struct_array(DDR_list_struct.(1))
DDR_file_struct_ident = mro_crism_parse_identifier(filename_struct = DDR_file_struct)
DDR_indx = cmset_op(DDR_file_struct_ident, 'AND', C_list_struct.(1), /index)
DDR_input_file_list = DDR_list_struct.(1)[DDR_indx]
ddr_file_list = DDR_input_file_list

;PHT_list_file = '/disks/dc026/msp_mos/north_polar_mapping/vnir/PHT/.txt'
PHT_list_file = root + dtype + '/PHT/' + obs_mode + '/planum_boreum_' + obs_mode + '_' + dtype_short + '_TRR_PRE_PHT.txt'
PHT_list_struct = mro_crism_read_list(PHT_list_file)
PHT_file_struct = mro_crism_filename_struct_array(PHT_list_struct.(1))
PHT_file_struct_ident = mro_crism_parse_identifier(filename_struct = PHT_file_struct)
PHT_indx = cmset_op(PHT_file_struct_ident, 'AND', C_list_struct.(1), /index)
PHT_input_file_list = PHT_list_struct.(1)[PHT_indx]
trdr_file_list = PHT_input_file_list

;stop

if (keyword_set(restore_initial)) then begin
    print, 'Restoring: ' + initial_save
    restore, initial_save
    print, 'Done'
endif else begin

if ((n_elements(ddr_list_file) GT 0) AND (n_elements(trdr_list_file) GT 0)) then begin
    ddr_list_struct = mro_crism_read_list(ddr_list_file)
    ddr_file_list = ddr_list_struct.active_array
    num_ddr_file_list = ddr_list_struct.num_active_array

    trdr_list_struct = mro_crism_read_list(trdr_list_file)
    trdr_file_list = trdr_list_struct.active_array
    num_trdr_file_list = trdr_list_struct.num_active_array

    ;stop

endif

help, ddr_file_list
help, trdr_file_list

ddr_file_struct_array = mro_crism_filename_struct_array(ddr_file_list)
help, ddr_file_struct_array
help, ddr_file_struct_array, /str

trdr_file_struct_array = mro_crism_filename_struct_array(trdr_file_list)
help, trdr_file_struct_array
help, trdr_file_struct_array, /str

ddr_identifier = mro_crism_parse_identifier(filename_struct = ddr_file_struct_array)
trdr_identifier = mro_crism_parse_identifier(filename_struct = trdr_file_struct_array)

trdr_samples_vec = make_array(n_elements(trdr_identifier), value = long(0))
trdr_lines_vec = make_array(n_elements(trdr_identifier), value = long(0))
trdr_bins_vec = make_array(n_elements(trdr_identifier), value = long(0))

ddr_match_indx = make_array(n_elements(trdr_identifier), value = long(-1))


;initial setup - sort out the scale of the problem

for t = 0, n_elements(trdr_identifier) -1 do begin
    trdr_indx = t
    ddr_indx = where(ddr_identifier EQ trdr_identifier[trdr_indx], num_ddr_indx)
    if (num_ddr_indx NE 1) then begin
      print, 'something wicked this way comes...'
  stop
    endif
    ddr_match_indx[t] = ddr_indx
   
    trdr_label_struct = mro_crism_quick_read(trdr_file_struct_array[trdr_indx].filename, /label_struct)
    ddr_label_struct = mro_crism_quick_read(ddr_file_struct_array[ddr_indx].filename, /label_struct)

    help, trdr_label_struct, /struct
    help, ddr_label_struct, /struct

    trdr_detector_mask = mro_crism_detector_mask(trdr_label_struct, /cross_track)
    trdr_detector_mask_total = total(trdr_detector_mask)

    trdr_samples_vec[t] = trdr_detector_mask_total
    trdr_lines_vec[t] = trdr_label_struct.lines
    
    trdr_bins_vec[t] = ceil(float(trdr_lines_vec[t]) / trdr_samples_vec[t])

    if (t EQ 0) then begin
      num_bands = long(trdr_label_struct.bands)
      lat_indx = where(strpos(strupcase(ddr_label_struct.band_name), 'LATITUDE, AREOCENTRIC') NE -1)
      lon_indx = where(strpos(strupcase(ddr_label_struct.band_name), 'LONGITUDE, AREOCENTRIC') NE -1)
  
  inc_indx = where(strpos(strupcase(ddr_label_struct.band_name), 'INA AT AREOID') NE -1)
  emi_indx = where(strpos(strupcase(ddr_label_struct.band_name), 'EMA AT AREOID') NE -1)
  pha_indx = where(strpos(strupcase(ddr_label_struct.band_name), 'PHASE ANGLE') NE -1)
  
  det_mask_min_indx = min(where(trdr_detector_mask EQ 1))
  det_mask_max_indx = max(where(trdr_detector_mask EQ 1))
    endif
    ;stop
endfor  ;t

print, 'Saving: ' + initial_save
save, /variables, filename = initial_save

endelse ;initial

;stop

vector_save = save_path + 'proximal_vectors_' + lower_longi + '_' + upper_longi + '.sav'

if (keyword_set(restore_vectors)) then begin
    print, 'Restoring: ' + vector_save
    restore, vector_save
    print, 'Done'
endif else begin

;init bin vectors
bin_total = total(trdr_bins_vec)

bin_count = make_array(1, bin_total, value = long(-1))
bin_identifier = make_array(2, bin_total, value = long(-1))

bin_latitude = make_array(1, bin_total, value = float(0.0))
bin_longitude = make_array(1, bin_total, value = float(0.0))

bin_incidence = make_array(1, bin_total, value = float(0.0))
bin_emission = make_array(1, bin_total, value = float(0.0))
bin_phase = make_array(1, bin_total, value = float(0.0))

bin_mean = make_array(num_bands, bin_total, value = float(-1.0))
bin_stddev = make_array(num_bands, bin_total, value = float(-1.0))

bin_flag = make_array(num_bands, bin_total, value = byte(0))

;populate bin vectors

mid_strip_distance = make_array(n_elements(trdr_identifier), value = float(0.0))

bin_runner = long(0)

for t = 0, n_elements(trdr_identifier) -1 do begin
    print
    print, '*** ' + strtrim(string(t+1),2) + ' / ' + strtrim(string(n_elements(trdr_identifier)),2) + ' ***'

    trdr_indx = t
    trdr_struct = mro_crism_quick_read(trdr_file_struct_array[trdr_indx].filename, /load_struct)
    ddr_struct = mro_crism_quick_read(ddr_file_struct_array[ddr_match_indx[trdr_indx]].filename, /load_struct)

    ;stop

    flag_total = total(trdr_struct.image EQ -0.10)
    if (flag_total GT 0) then begin
      print, '$$$'
      print, trdr_file_struct_array[trdr_indx].file_name
  print, flag_total
  print, '$$$'
    endif
    
    trdr_size = size(trdr_struct.image, /struct)
    ddr_size = size(ddr_struct.image, /struct)

    help, trdr_struct, /struct
    help, ddr_struct, /struct

    mid_strip_dist = great_circle_distance(ddr_struct.image[det_mask_min_indx, ddr_size.dimensions[1]/2, lat_indx], ddr_struct.image[det_mask_min_indx, ddr_size.dimensions[1]/2, lon_indx], $
                                 ddr_struct.image[det_mask_max_indx, ddr_size.dimensions[1]/2, lat_indx], ddr_struct.image[det_mask_max_indx, ddr_size.dimensions[1]/2, lon_indx], $
             /double)

    mid_strip_distance[t] = mid_strip_dist
    
    for s = 0, trdr_bins_vec[t] -1 do begin
      ;lat_sample = ddr_struct.image[det_mask_min_indx:det_mask_max_indx, (s * trdr_samples_vec[s]):((s + 1) * trdr_samples_vec[s] -1) < (ddr_size.dimensions[1] -1), lat_indx]
      ;lon_sample = ddr_struct.image[det_mask_min_indx:det_mask_max_indx, (s * trdr_samples_vec[s]):((s + 1) * trdr_samples_vec[s] -1) < (ddr_size.dimensions[1] -1), lon_indx]
      lat_sample = ddr_struct.image[det_mask_min_indx:det_mask_max_indx, (s * trdr_samples_vec[t]):((s + 1) * trdr_samples_vec[t] -1) < (ddr_size.dimensions[1] -1), lat_indx]
      lon_sample = ddr_struct.image[det_mask_min_indx:det_mask_max_indx, (s * trdr_samples_vec[t]):((s + 1) * trdr_samples_vec[t] -1) < (ddr_size.dimensions[1] -1), lon_indx]

      lat_sample_size = size(lat_sample, /struct)
      lon_sample_size = size(lon_sample, /struct)

      lat_sample_ref = lat_sample[lat_sample_size.dimensions[0]/2, lat_sample_size.dimensions[1]/2]
      lon_sample_ref = lon_sample[lon_sample_size.dimensions[0]/2, lon_sample_size.dimensions[1]/2]

      bin_latitude[bin_runner] = lat_sample_ref
      bin_longitude[bin_runner] = lon_sample_ref

      bin_count[bin_runner] = lat_sample_size.n_elements
      bin_identifier[*,bin_runner] = [t,s]

      ;ina_sample = ddr_struct.image[det_mask_min_indx:det_mask_max_indx, (s * trdr_samples_vec[s]):((s + 1) * trdr_samples_vec[s] -1) < (ddr_size.dimensions[1] -1), inc_indx]
      ;emi_sample = ddr_struct.image[det_mask_min_indx:det_mask_max_indx, (s * trdr_samples_vec[s]):((s + 1) * trdr_samples_vec[s] -1) < (ddr_size.dimensions[1] -1), emi_indx]
      ;pha_sample = ddr_struct.image[det_mask_min_indx:det_mask_max_indx, (s * trdr_samples_vec[s]):((s + 1) * trdr_samples_vec[s] -1) < (ddr_size.dimensions[1] -1), pha_indx] 
      ina_sample = ddr_struct.image[det_mask_min_indx:det_mask_max_indx, (s * trdr_samples_vec[t]):((s + 1) * trdr_samples_vec[t] -1) < (ddr_size.dimensions[1] -1), inc_indx]
      emi_sample = ddr_struct.image[det_mask_min_indx:det_mask_max_indx, (s * trdr_samples_vec[t]):((s + 1) * trdr_samples_vec[t] -1) < (ddr_size.dimensions[1] -1), emi_indx]
      pha_sample = ddr_struct.image[det_mask_min_indx:det_mask_max_indx, (s * trdr_samples_vec[t]):((s + 1) * trdr_samples_vec[t] -1) < (ddr_size.dimensions[1] -1), pha_indx]  


      ina_sample_ref = ina_sample[lat_sample_size.dimensions[0]/2, lat_sample_size.dimensions[1]/2]
      emi_sample_ref = emi_sample[lat_sample_size.dimensions[0]/2, lat_sample_size.dimensions[1]/2]
      pha_sample_ref = pha_sample[lat_sample_size.dimensions[0]/2, lat_sample_size.dimensions[1]/2]

      bin_incidence[bin_runner] = ina_sample_ref
      bin_emission[bin_runner] = emi_sample_ref
      bin_phase[bin_runner] = pha_sample_ref
    
      ;trdr_sample = trdr_struct.image[det_mask_min_indx:det_mask_max_indx, (s * trdr_samples_vec[s]):((s + 1) * trdr_samples_vec[s] -1) < (ddr_size.dimensions[1] -1), *]
      trdr_sample = trdr_struct.image[det_mask_min_indx:det_mask_max_indx, (s * trdr_samples_vec[t]):((s + 1) * trdr_samples_vec[t] -1) < (ddr_size.dimensions[1] -1), *]
      trdr_sample_stats = mro_crism_image_stats(trdr_sample, /simple, /sigma)

      ;stop
      
   bin_mean[*,bin_runner] = trdr_sample_stats.image_mean
   bin_stddev[*,bin_runner] = trdr_sample_stats.image_stddev

      ;bin_flag[*,bin_runner] = (total(total(trdr_sample EQ -0.10, 1), 1) GT 0)
  
      ;if (total((finite(trdr_sample_stats.image_mean[52]) EQ 0)) GT 0) then begin
  ;    ;stop
      ;endif
  
  ;if (total((finite(trdr_sample_stats.image_stddev[52]) EQ 0)) GT 0) then begin
  ;    ;stop
  ;endif
  
      ;stop

      bin_runner++
    endfor

    ;stop

endfor  ;t

print, 'Saving: ' + vector_save
save, mid_strip_distance, bin_total, bin_count, bin_identifier, bin_latitude, bin_longitude, bin_incidence, bin_emission, bin_phase, bin_mean, bin_stddev, bin_flag, filename = vector_save

endelse ;vector

matrix_save = save_path + 'weighted_proximal_matricies_distance_and_weight_' + lower_longi + '_' + upper_longi +  '_' + band + '.sav'
matrix_save_specific_band = save_path + 'weighted_proximal_matricies_' + lower_longi + '_' + upper_longi + '_' + band + '.sav'

if (keyword_set(restore_matricies)) then begin
    print, 'Restoring: ' + matrix_save
    restore, matrix_save
    restore, matrix_save_specific_band
endif else begin

;establish distance and count relationship matricies
;establish mean-delta and stddev-ratio matricies

;custom_indx = [1, 2, 3, 8]
;custom_indx = [1]; 440
;custom_indx = [2]; 530
;custom_indx = [3]; 600
;custom_indx = [8]; 770

bands = transpose(trdr_label_struct.band_name[custom_indx])
;specific_bands = individual_bands[1:3] ;still need to incorporate band #8
num_bands_matrix = n_elements(bands)

bin_latitude_radians = bin_latitude * !dtor
bin_longitude_radians = bin_longitude* !dtor

distance_count_matrix = make_array(bin_total, bin_total, value = -1.0, /float)
stats_compare_matrix = make_array(bin_total, bin_total, num_bands_matrix, value = !values.f_nan, /float)

matrix_sample = lindgen(bin_total)
matrix_line = lindgen(bin_total)
matrix_sample_matrix = cmreplicate(matrix_sample, bin_total)
matrix_line_matrix = transpose(cmreplicate(matrix_line, bin_total))

ld_indx = where(matrix_line_matrix GT matrix_sample_matrix, num_ud_indx)
ud_indx = where(matrix_line_matrix LT matrix_sample_matrix, num_ld_indx)

for i = 0, bin_total -1 do begin
    ;print
    ;print, '### ' + strtrim(string(i+1),2) + ' / ' + strtrim(string(bin_total),2) + ' ###'
    
    for j = i + 1, bin_total -1 do begin

      sample_distance = great_circle_distance(bin_latitude_radians[i], bin_longitude_radians[i], bin_latitude_radians[j], bin_longitude_radians[j], /radians, /double)
      distance_count_matrix[i,j] = sample_distance
  
      sample_count_geometric_mean = (bin_count[i] * bin_count[j])^(1.0/2.0)
      distance_count_matrix[j,i] = sample_count_geometric_mean
      
      ;mean_delta = bin_mean[*,i] - bin_mean[*,j]
      ;stddev_ratio = bin_stddev[*,i] / bin_stddev[*,j]
      ;stats_compare_matrix[i,j,*] = mean_delta
      ;stats_compare_matrix[j,i,*] = stddev_ratio
      mean_delta = bin_mean[custom_indx,i] - bin_mean[custom_indx,j]
      stddev_ratio = bin_stddev[custom_indx,i] / bin_stddev[custom_indx,j]
      stats_compare_matrix[i,j,custom_indx] = mean_delta
      stats_compare_matrix[j,i,custom_indx] = stddev_ratio
      
    
    ;  if (finite(mean_delta[52]) EQ 0) then begin
          ;stop
    ;  endif
  
 ; if (finite(stddev_ratio[52]) EQ 0) then begin
      ;stop
 ; endif
    
   endfor  ;j
endfor  ;i

;mro_crism_quick_view, distance_count_matrix, winid = 10, min = 0, max = max(distance_count_matrix[ld_indx]), resize= 0.250, /flip
;mro_crism_quick_view, distance_count_matrix, winid = 11, min = min(distance_count_matrix[ud_indx]), max = max(distance_count_matrix[ud_indx]), resize= 0.250, /flip

;establish distance gaussian weighting
weighting_sigma = 2.0 * median(mid_strip_distance)
;weighting_sigma = 4.0 * median(mid_strip_distance)

weighting_abscissa = findgen(5.0e5 + 1) / 1.0e6
weighting_function = scaled_offset_normal_distribution(weighting_abscissa, [0.0, weighting_sigma, 1.0, 0.0])
weighting_function = weighting_function / max(weighting_function)

weighting_function_spline_init = spl_init(weighting_abscissa, weighting_function)

distance_abscissa = distance_count_matrix[ld_indx]
distance_weight = spl_interp(weighting_abscissa, weighting_function, weighting_function_spline_init, distance_abscissa)

weighting_matrix = distance_count_matrix * 0.0
;mro_crism_quick_view, weighting_matrix, winid = 31, min = 0, max = 1, resize= 0.250, /flip
;stop

weighting_matrix[ld_indx] = distance_weight
;mro_crism_quick_view, weighting_matrix, winid = 12, min = 0, max = 1, resize= 0.250, /flip

weighting_matrix[ud_indx] = (transpose(weighting_matrix) * distance_count_matrix)[ud_indx]
;mro_crism_quick_view, weighting_matrix, winid = 13, min = min(weighting_matrix[ud_indx]), max = max(weighting_matrix[ud_indx]), resize= 0.250, /flip

;mro_crism_quick_view, (weighting_matrix GT 0), winid = 14, min = 0, max = 1, resize= 0.250, /flip

stop

print, 'Saving: ' + matrix_save
save, distance_count_matrix, weighting_matrix, ld_indx, ud_indx, filename = matrix_save
print, 'Saving: ' + matrix_save_specific_band
save, stats_compare_matrix, filename = matrix_save_specific_band
;save, distance_count_matrix, weighting_matrix, ld_indx, ud_indx, filename = matrix_save

endelse ;matrix


stop

  mro_crism_proximal_optimization_linear_system_lim, distance_count_matrix = distance_count_matrix, bin_identifier = bin_identifier, stats_compare_matrix = stats_compare_matrix, $ ;weighting_matrix, 
                        bin_identifier = bin_identifier, $ bin_flag ;stats_compare_matrix, $
  hold_indx = hold_indx


END
                                                                                                                                                                                                                                                                                                                                                                                                                                                   parse_my_info_lim.pro                                                                               0000664 €    80002301 00000010026 13144361573 014574  0                                                                                                    ustar   limh1                           crism-soc                                                                                                                                                                                                              ;07/14/2017 (limh1)
;    -Parse out attribute information by MY and Ls
;

PRO parse_my_info, dtype = dtype, obs_mode = obs_mode

Working_directory = '/disks/dc026/msp_mos/north_polar_mapping/' + dtype + '/SLF/PRJ/' + obs_mode + '_attributes'
print, 'Working directory: ' + working_directory

;==== Potentially relevant files ====;
;attribute_files = file_search('/disks/dc026/msp_mos/north_polar_mapping/vnir/SLF/PRJ/HSP_attributes/*.txt', count = num_attr_files)
;attributes_data_array = PRJ_directory + 'attributes_data_array.csv'

A1_list_file = working_directory + '/filenames.txt'
A2_list_file = working_directory + '/classtypes.txt'
A3_list_file = working_directory + '/targetIDs.txt'
A4_list_file = working_directory + '/segmentIDs.txt'
A5_list_file = working_directory + '/starttimes.txt'
A6_list_file = working_directory + '/slongis.txt'
A7_list_file = working_directory + '/detector.txt'

A5_MY_list_file = working_directory + '/starttimes_MY_really.txt' ;A5 in MY

;==== READ COLUMN OF ATTRIBUTES ====;
readcol, A1_list_file, FullFilename, format = 'A'
readcol, A2_list_file, Classtype, format = 'A'
readcol, A3_list_file, targetID, format = 'A'
readcol, A4_list_file, segmentID, format = 'A'
readcol, A5_MY_list_file, MY, format = 'I'
readcol, A6_list_file, SolarLongi, format = 'F'
readcol, A7_list_file, A7list, format = 'A'

;==== Num check if desired ====;
num_ctypes = n_elements(Classtype)
num_tIDs = n_elements(targetID)
num_sIDs = n_elements(segmentID)
num_MYs = n_elements(MY)
num_Ls = n_elements(SolarLongi)

  min_my = min(MY)
  max_my = max(MY)

  ls_bounds = lindgen(37) * 10

;stop

;==== LOG INFO INTO TEXTFILES ====; 
  for i = 1, n_elements(ls_bounds) -1 do begin
    
    Ls_window_list = ''
    
     for j = min_my, max_my do begin
     
        index = where((MY EQ j) and (solarlongi GE ls_bounds[i-1]) AND (SolarLongi LT ls_bounds[i]), num_index)
      
       if num_index GT 0 then begin
        print, j, ls_bounds[i], ls_bounds[i-1]
        print, num_index

        print, Classtype[index] + targetID[index] + '_' + segmentID[index]
        working_file = working_directory + '/MY/' + '/MY' + strtrim(string(j), 2) + '_Ls_' + strtrim(string(ls_bounds[i-1], format = '(I3.3)'), 2) + '_' + strtrim(string(ls_bounds[i], format = '(I3.3)'), 2) + '.txt'
        print, working_file

        print, FullFilename[index] ;For ddr list
        working_file_full_DDR = working_directory + '/MY/' + '/MY' + strtrim(string(j), 2) + '_Ls_' + strtrim(string(ls_bounds[i-1], format = '(I3.3)'), 2) + '_' + strtrim(string(ls_bounds[i], format = '(I3.3)'), 2) + '_DDR_paths.txt'
        print, working_file_full_DDR

;stop
      
        openw, out_lun_all, working_file, /get_lun
        ;printf, out_lun_all, j, ls_bounds[i-1], ls_bounds[i]
        printf, out_lun_all, num_index
        printf, out_lun_all, transpose(Classtype[index] + targetID[index] + '_' + segmentID[index])
        free_lun, out_lun_all
        
        openw, out_lun_ddr, working_file_full_DDR, /get_lun
        printf, out_lun_ddr, num_index
        printf, out_lun_ddr, transpose(FullFilename[index])
        free_lun, out_lun_ddr

          if (Ls_window_list[0] EQ '') then begin
            Ls_window_list = [Classtype[index] + targetID[index] + '_' + segmentID[index]]
            ;Ls_window_list = [FullFilename[index]]
          endif else begin
              Ls_window_list = [Ls_window_list, [Classtype[index] + targetID[index] + '_' + segmentID[index]]]         
              ;Ls_window_list = [Ls_window_list, [FullFilename[index]]]
            endelse
      
       endif
     
     
     endfor



    if (strlen(Ls_window_list[0]) GT 0) then begin
     other_working_file = working_directory + '/LS/' + 'Ls_' + strtrim(string(ls_bounds[i-1], format = '(I3.3)'), 2) + '_' + strtrim(string(ls_bounds[i], format = '(I3.3)'), 2) + '.txt'
     print, other_working_file

       openw, out_lun_Ls, other_working_file, /get_lun
       printf, out_lun_Ls, n_elements(Ls_window_list)
       printf, out_lun_Ls, transpose(Ls_window_list)
       free_lun, out_lun_Ls
       
    endif
  
  
  endfor


END                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          write_mosaic_file_lim.pro                                                                           0000664 €    80002301 00000003062 13144343337 015426  0                                                                                                    ustar   limh1                           crism-soc                                                                                                                                                                                                              ;06/27/17 (limh1)
;   -Writing out mosaic files for mapping
;

pro write_mosaic_file_lim, output_file = output_file, input_files = input_list_file, input_fids = input_fids, $
  input_dims_array = input_dims_array, input_bnames = input_bnames, mosaic_xsize = mosaic_xsize, $
  mosaic_ysize = mosaic_ysize, input_pxl_size = input_pxl_size, mosaic_x0 = mosaic_x0, mosaic_y0 = mosaic_y0
  
xsize_pix = round(mosaic_xsize/input_pxl_size[0])
ysize_pix = round(mosaic_ysize/input_pxl_size[1])

  num_file_line_array = 5 * n_elements(input_list_file) + 3.
  header_line_array = ['ENVI MOSAIC TEMPLATE (G)', 'Output Size: ' + strtrim(string(xsize_pix), 1) + ' x ' + strtrim(string(ysize_pix), 1)]
  openw, out_mos_lun, output_mosaic_file, /get_lun
  
  for i = 0,  n_elements(input_list_file) -1 do begin
    
    if (i EQ 0) then begin
    
      for k = 0, n_elements(header_line_array) -1 do begin
        print, header_line_array[k]
        printf, out_mos_lun, header_line_array[k]
    
      endfor
    
      print
      printf, out_mos_lun
    
    endif
    
   ; printf, out_mos_lun, 'File : '+input_file_list[i]
    printf, out_mos_lun, 'File : '+input_list_file[i]
    printf, out_mos_lun, 'Bands: 1-'+strtrim(string(n_elements(input_bnames)),1)
    printf, out_mos_lun, 'Dims : 1-'+strtrim(string(input_dims_array[2,i]+1),1)+',1-'+ strtrim(string(input_dims_array[4,i]+1),1)
    printf, out_mos_lun, 'Info : ('+strtrim(string(mosaic_x0[i]), 1) +',' + strtrim(string(mosaic_y0[i]), 1)+') {0} [ {0}] {} ^1^2.0^'
    printf, out_mos_lun
    
  endfor
  
  free_lun, out_mos_lun
  
end
                                                                                                                                                                                                                                                                                                                                                                                                                                                                              write_my_attributes_lim.pro                                                                         0000664 €    80002301 00000003217 13144350573 016051  0                                                                                                    ustar   limh1                           crism-soc                                                                                                                                                                                                              ;07/30/2017 (limh1)
;    -Write attributes from shapefiles into their own textfiles
;-


PRO write_my_attributes, dtype = dtype, obs_mode = obs_mode, A1 = A1, A2 = A2, A3 = A3, A4 = A4, A5 = A5, A6 = A6, A7 = A7, $
                         filename = filename, classtype = classtype, targetID = targetID, segmentID = segmentID, $
                         starttime = starttime, slongi = slongi, detector = detector


PRJ_directory = '/disks/dc026/msp_mos/north_polar_mapping/' + dtype + '/SLF/PRJ/' + obs_mode + '_attributes/'

;A1
A1_list_file = PRJ_directory + 'filenames.txt'
openw, unit1, A1_list_file, /get_lun
printf, unit1, n_elements(A1)
printf, unit1, transpose(A1)
free_lun, unit1

;A2
A2_list_file = PRJ_directory + 'classtypes.txt'
openw, unit2, A2_list_file, /get_lun
printf, unit2, n_elements(A2)
printf, unit2, transpose(A2)
free_lun, unit2

;A3
A3_list_file = PRJ_directory + 'targetIDs.txt'
openw, unit3, A3_list_file, /get_lun
printf, unit3, n_elements(A3)
printf, unit3, transpose(A3)
free_lun, unit3

;A4
A4_list_file = PRJ_directory + 'segmentIDs.txt'
openw, unit4, A4_list_file, /get_lun
printf, unit4, n_elements(A4)
printf, unit4, transpose(A4)
free_lun, unit4

;A5
A5_list_file = PRJ_directory + 'starttimes.txt'
openw, unit5, A5_list_file, /get_lun
printf, unit5, n_elements(A5)
printf, unit5, transpose(A5)
free_lun, unit5

;A6
A6_list_file = PRJ_directory + 'slongis.txt'
openw, unit6, A6_list_file, /get_lun
printf, unit6, n_elements(A6)
printf, unit6, transpose(A6)
free_lun, unit6

;A7
A7_list_file = PRJ_directory + 'detector.txt'
openw, unit7, A7_list_file, /get_lun
printf, unit7, n_elements(A7)
printf, unit7, transpose(A7)
free_lun, unit7

END                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
