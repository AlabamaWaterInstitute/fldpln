#  NOTIFICATION
Data is available in an archive on Hydroshare. :arrow_right:[Hydroshare](https://www.hydroshare.org/resource/2ba43947ef6447beaf055349c883c96e/)<br>
Data in this folder will support an example of producing an FLDPLN library.<br>
<br>
<br>
<br>

#  INFORMATION
##  DEVICE
###  WORKSTATION
*  Processor:<br>
    Intel(R) Xeon(R) CPU E5-2630 (2 processors)<br>
*  Installed RAM:<br>
    128 GB DDR-3 1333MHz (10 of 16)
*  System type:<br>
    64-bit operating system, x64-based processor
*  Edition:<br>
    Windows Server 2022 Standard
*  Version:<br>
    21H2
*  OS Build:<br>
    20348.1787
###  LAPTOP
*  Processors:<br>
    12th Gen Intel(R) Core(TM) i7-12700H 2.3 GHz
*  Installed RAM:<br>
    32 GB
*  System type:<br>
    64-bit operating system, x64-based processor
*  Edition:<br>
    Windows 10 Education
*  Version:<br>
    21H2
*  OS build:<br>
    19044.3086
*  Experience:<br>
    Windows Feature Experience Pack

##  MATLAB
*  Version:<br>
    R2022b Update 5
*  Type:<br>
    64-bit (win64)
*  Parallel Tools:<br>
    Parallel Computation Tools from Matlab ADD-Ons

##  ESRI GIS
###  ArcGIS
*  Version:<br>
    ArcGIS Pro 3.1.2
###  ArcHydro Tool		
*  Version:<br>
    ArcGIS Pro 3.1 :arrow_right:[HydroTool](https://www.esri.com/en-us/industries/water-resources/arc-hydro/downloads)

##  DATA
###  HYDRO UNIT CODE (HUC)
*  HUC-8 :arrow_right:[HUC](https://www.hydroshare.org/resource/b832a6c2f96541808444ec9562c5247e/) or :arrow_right:[WDB](https://hydro.nationalmap.gov/arcgis/rest/services/wbd/MapServer)
###  DIGITAL ELEVATION MODEL (DEM)
*  10-meter Resolution :arrow_right:[USGS National Map](https://apps.nationalmap.gov/downloader/)

##  MODEL
###  Valley Floor Mapper (VFM)
*  Request to Dr. Jude Kastens at the University of Kansas
###  FLDPLN Model
*  Request to Dr. Jude Kastens at the University of Kansas
<br>
<br>
<br>

#  1. PREPROCESSING
*  The Amite River Basin in Louisiana (LA.) is an example of preprocessing.<br>
*  To initiate the study, download the HUC-8 Watershed Boundary Data (WBD) and the DEM of the study area.<br>
*  The HUC-8 refers to a shapefile that outlines the watershed boundary of the United States (U.S.).<br>
*  The DEM refers to the elevation data of the United States, sourced from the United States Geological Survey (USGS).<br>
*  To obtain DEM data from the USGS, users can upload their own shapefile to the website.<br>
*  However, keep in mind that there may be unexpected errors, as mentioned on the website.<br>

##  [1]  EXTRACTING A SHAPEFILE OF THE WATERSHED FOR THE STUDY AREA
To distinguish the designated study area, the HUC-8, which is the Continental United States (CONUS) WBD as the type of polygon, is necessary.<br>

    [Step 1]  Open the HUC-8 shapefile in ArcGIS Pro.
    [Step 2]  Right-click the shapefile and go to the attribute table.
    [Step 3]  Delete data in the attribute table that does not relate to the target study area.
    [Step 4]  Right-click the shapefile, click the Data, and choose Export Features to save.
              Save Ex. Amite_basin.shp
    [Result]  You will only be provided with the shapefile for the Amite River watershed boundary.

![F1_HUC_8_polygon](https://github.com/Jsong33AWI/FLDPLN/assets/122640336/3495262c-584c-41e8-b86c-6d436885755a)

![F2_Extracted_AmiteRiver](https://github.com/Jsong33AWI/FLDPLN/assets/122640336/0d958d44-0fb2-451a-b55f-34f15a181cb5)

##  [2]  COMBINING RAW DEMs
If the study area covers more than one piece of DEM, it is necessary to combine them through processing.<br>
When conducting a large-scale study of FIM, it is common to combine DEMs.<br>
However, this process is not necessary if the study area is already located on a single DEM.<br>

    [Step 1]  In the Data Management Tool in ArcGIS Pro, go to the Mosaic to New Raster.
    [Step 2]  Type the information into Input Rasters, Output Location, Extension, and Number of Bands. 
              The other options are not considered.
              The Amite River Basin covers 4 pieces of raw DEM, thus the 4 DEMs should be placed into the Input Rasters.
              Save Ex. Combined_DEM.tif 
    [Result]  Once the 4 pieces of the downloaded DEM are combined, you will have a single, unified DEM.

![F3_CombinedDEM](https://github.com/Jsong33AWI/FLDPLN/assets/122640336/da5eb567-5c2d-4b19-8be6-c4112cf151fc)

![F4_UnifiedDEM](https://github.com/Jsong33AWI/FLDPLN/assets/122640336/6d216f8c-5fa4-4ad9-a511-e8439816c0ab)

##  [3]  MASKING THE COMBINED DEM [2] USING EXTRACTED SHAPEFILE [1]
By applying the created shapefile to mask the DEM, users can obtain a masked DEM for the study area only.<br>
The square box DEM, which is a raw elevation model from the USGS, contains data on elevations that include unnecessary locations.<br>
The following steps are showing how to mask the study area DEM by cutting the unnecessary part of the DEM off.<br>

    [Step 1]  Open the extracted shapefile (Amite_basin.shp) from [1] and the DEM (Combined_DEM.tif) from [2].
    [Step 2]  At the Spatial Analysis Tools, go to the Extract by Mask.
    [Step 3]  The DEM and extracted shapefile should be placed into the input raster and feature mask, respectively. 
    [Step 4]  Make sure to include the file extension as TIF when naming the output file.
              Save Ex. Masked_Amite_DEM.tif
    [Result]  It appears that the masked DEM and extracted shapefile are matched.

![F5_MaskingProcess_Masking](https://github.com/Jsong33AWI/FLDPLN/assets/122640336/9bcf9f7e-972c-40d3-b4af-ac2b96ba1064)

##  [4]  FILL DEM
    [Step 1]  Go to the Terrain Preprocessing in ArcHydro Tool.
    [Step 2]  In the DEM Manipulation, click the Fill Sinks.
    [Step 3]  In the Fill Sinks, Type the name and save it in TIF format.
              Save Ex. Fil.tif
    [Result]  You will see the Filled DEM map in the main GIS window.

![F6_FilDEM](https://github.com/Jsong33AWI/FLDPLN/assets/122640336/5681e2af-e4dc-40b4-892d-7ff6598e9312) 

##  [5]  FLOW DIRECTION
    [Step 1]  Go to the Terrain Preprocessing in ArcHydro Tool.
    [Step 2]  Click the Flow Direction.
    [Step 3]  Input is the output from Fill DEM (Fil.tif).
    [Step 4]  Type the name and save it in TIF format.
              Save Ex. Fdr.tif
    [Result]  You will see the Flow Direction map in the main GIS window.

![F7_Fdr](https://github.com/Jsong33AWI/FLDPLN/assets/122640336/5400e345-ecc2-4f90-8a9a-718522d2accb)

##  [6]  FLOW ACCUMULATION
    [Step 1]  Go to the Terrain Preprocessing in ArcHydro Tool.
    [Step 2]  Click the Flow Accumulation.
    [Step 3]  Input is the output from Flow Direction (Fdr.tif).
    [Step 4]  Type the name and save it in TIF format.
              Save Ex. Fac.tif
    [Result]  You will see the Flow Accumulation map in the main GIS window.

##  [7]  CONVERTING FORMAT FROM TIF TO ESRI BIL
    [Step 1]  Go to the ArcGIS pro, Catalog pane.
    [Step 2]  In the Catalog pane, click the Computer tab and go to the file location you saved.
              Ex. Masked_Amite_DEM.tif, Fil.tif, Fdr.tif, Fac.tif
    [Step 3]  Right-click on each TIF file and click the export to a different format.
    [Step 4]  Convert all TIF files to ESRI BIL format. 
    [Step 5]  During the conversion of files, other options are not considered.
              Only for Fac.tif converting, type the -9999 in NoData Value.

![F8_Fac](https://github.com/Jsong33AWI/FLDPLN/assets/122640336/275855f7-ecca-4752-b7b6-67d34ba5e6ae)

##  [8]  DATA AFTER THE PROCESSING
After completing this preprocessing, you should have 8 items listed below.<br>
1. Masked_Amite_DEM.tif and Masked_Amite_DEM.bil<br>
2. Fil.tif and Fil.bil<br>
3. Fdr.tif and Fdr.bil<br>
4. Fac.tif and Fac.bil<br>
<br>
<br>
<be>

#  2. FLDPLN MODEL
##  [1]  floodplain_segdb_Tuscaloosademo.m
*  Home directory: 15 - 20 line

        >    dr1 = 'Tuscaloosa\';
        >    dr0 = ['C:\Users\jsong33\Desktop\Work\FIM\FLDPLN\StudyAreas\Tuscaloosa\',dr1];
        >    if(~exist(dr0,'dir'))
        >    dr0 = ['C:\Users\jsong33\Desktop\Work\FIM\FLDPLN\StudyAreas\Tuscaloosa\',dr1];
        >    end

*  Iteration step size (in vertical DEM units): 27 - 32 line

        >    %dh = 0.01; 
        >    %dh = 0.1;
        >    %dh = 0.25;
        >    %dh = 0.5;
        >    %dh = 1;

*  Identify synthetic stream network: 54 - 58 line

        >    mi2px = 5280*.3048/pxsz; % horizontal units assumed to be meters (otherwise modify as needed) 
        >    seglen = round(5*mi2px); % 5 miles maximum segment length
        >    sqmi2px = (5280^2)*(.3048^2)/(pxsz^2); % horizontal units assumed to be meters (otherwise modify as needed)
        >    facthr = round(25*sqmi2px); % 25 sq. mi minimum catchment size (min for segment break)
        >    strthr = round(70*sqmi2px); % 70 sq. mi minimum catchment size (min for defining a stream)

*  Seedpts: 69 - 70 line

        >    create_segdb_facthr_maxlen(strthr,facthr,seglen,fdrf,facf,segf,matf)          % With no seedpts
        >    create_segdb_facthr_maxlen_seedpts(shpf,facthr,seglen,fdrf,facf,segf,matf)    % With seedpts

*  Parallelization: 187 - 212 line 

        *  Worker Pool

        >    if(1>0)
        >        g = gcp('nocreate');
        >        if(~isempty(g))
        >            p = gcp;
        >        else
        >            p  = parpool;
        >        end
        >        for j = 1:num
        >            pf(j) = parfeval(p,@fldpln_model_v5ram,1,seg_list(j),inp,fildat,fdrdat,fldmn,fldmx,dh,mxht,1,bg);
        >        end
        >    %   wait(pf)
        >        pfresults = cell(num,1);
        >        for idx = 1:num
        >    %   fetchNext blocks until next results are available.
        >            [completedIdx,value] = fetchNext(pf);
        >            pfresults{completedIdx} = value;
        >            fprintf('Segment %d completed\n', completedIdx);
        >        end
        >    %   delete(gcp('nocreate'))
        >    end
        >    return;

        *  Simple One-segment-at-a-time
           1. Active lines below
           2. Change 1>0 to 1<0 in the IF statement from Worker Pool
            
        >    %    parfor j = 1:num
        >    %    fldpln_model_v5ram(seg_list(j),inp,fildat,fdrdat,fldmn,fldmx,dh,mxht,1,bg);
        >    %    end
        >    %    return;

##  [2]  create_fldpln_map_cks_Tuscaloosademo.m



















