# FLDPLN
The FLDPLN concept was documented by [Jude Kastens' 2008 Ph.D. dissertation](https://kuscholarworks.ku.edu/handle/1808/5354) and was, in that document, put forward as
> "... a new, physically-based computational model (called FLDPLN, or "Floodplain") for mapping potential inundation extents (floodplains) using gridded topographic data"

And further,
> "Due to the parametric economy of FLDPLN, this model has significant advantages over existing methods such as hydrodynamic models."

Since that publication, FLDPLN has been used in a number of instances to estimate the inundated extent of floods in Kansas and elsewhere for both research and real-time emergency management purposes.

We present here the code to produce the basic FLDPLN libraries, along with some explanation of the preprocessing needed to prepare the inputs from the raw DEM. Links and information regarding how to use these libraries in the development of maps of inundated extent for specific events are forthcoming.


# Demonstration
Using FLDPLN consists of three general steps:<br>

1. Preparing a DEM and other basic input data,
2. Computing a FLDPLN library, and
3. Using the library for mapping.

The following instructions will focus on **step 2.** to allow you to use a pre-made set of DEM and input resources and create a simple test FLDPLN library.

We have prepared a simple demonstration case based on our favorite city -- Tuscaloosa, Alabama!<br>
The input files are available as a [Package on hydroshare.org](https://www.hydroshare.org/resource/2ba43947ef6447beaf055349c883c96e/) (free, requires registration) and the steps below outline the process of preparing a set of FLDPLN library files based on that package.

You will need either [Matlab](https://www.mathworks.com/products/matlab.html) or [Octave](https://octave.org/download) installed to execute these scripts. We are actively exploring ports of the code to additional languages. Currently, the parallel execution only works with Matlab. If you know how to make it work in Octave, please let us know!<br>
We have provided another [Readme](data/data_readme.md) explaining the preprocessing for the data, if you would like to follow those steps to prepare your own example data.

## Required software versions
### Matlab ...
  * Version:<br>
    R2022b Update 5
  * Parallel tool:<br>
    Parallel Computation Tools from Matlab ADD-Ons

### or OCTAVE
  * Version:<br>
    8.2.0

### ESRI ArcGIS
***(Only needed if you do not use our pre-packaged example)***
  * Version:<br>
    ArcGIS Pro 3.1.2
  * ArcHydro Tool:<br>
    ArcHydro Tools Pro 3.1


# Steps to Run Demo:
1. Install Matlab or Octave according the versions indicated above.
We used a Mac and found the hint here handy to configure an alias to run Matlab on the command line. Clone this repository to your computer.<br>
```git clone https://github.com/AlabamaWaterInstitute/fldpln```

2. Download the [Hydroshare resource](https://www.hydroshare.org/resource/2ba43947ef6447beaf055349c883c96e/) and unzip the `Tuscaloosa` folder into the same directory where you have cloned the code.

3. From the main repository directory, run the test:<br>
``` matlab -r "run('src/floodplain_segdb_Tuscaloosademo.m')" ```

# Notes:
## Testing
The primary script for executing the FLDPLN model is the `fldpln_model_v5ram` in the src folder. But we will use the wrapper script `floodplain_segdb_Tuscaloosademo.m` to set the input variables.

## FLDPLN Model Settings

**[1] fldpln_model_v5ram<br>**

  * Octave does not currently support the `-v7,3` file specification and you will need to change the versino to `-v7.0` to run in Octave.<br>
  Find these three lines to make the change:<br>

```Matlab
save([segdr,outf,'_tmp.mat'],'fldpln_info','fldpln_bdy','bdy','bdct','fldht','ct_tot','proctime','-v7.0');
...
save([segdr,outf,'.mat'],'fldpln','header','-v7.0');
...
save([segdr,outf,'.mat'],'fldpln','header','-v7.0');
```

<br>
<br>
<br>

**[2] floodplain_segdb_Tuscaloosademo.m<br>**

  * Home directory.<br>
This home directory can be adjusted. Right now, it is set for execution from the src folder, looking for the files in the adjacent Tuscaloosa folder, which is what you will get from the unzipped hydroshare package.<br>

```Matlab
dr1 = 'Tuscaloosa/';
dr0 = ['../',dr1];
if(~exist(dr0,'dir'))
dr0 = ['../',dr1];
end
```

  * The FLDPLN iteration step size will drastically change the nature of the library. The value is set in vertical DEM units. More documentation is coming about that. Choose an iteration step size.<br>

```Matlab
% dh = 0.01;
% dh = 0.1;
% dh = 0.25;
  dh = 0.5;
% dh = 1;
```

  * Maximum flood depth.<br>

```Matlab
fldmx = 3; (maximum flood depth)
```

  * Read input and write output.<br>

```Matlab
demf = [dr0,'bil/DEM.bil']; (read raw DEM)
filf = [dr0,'bil/Fil.bil']; (read filled DEM)
fdrf = [dr0,'bil/Fdr.bil']; (read flow direction)
facf = [dr0,'bil/Fac.bil']; (read flow accumulation)

strf0 = 'str.bil';
strf = [dr0,'bil/',strf0]; (write stream segments)
segf = [dr0,'bil/str_segid.bil']; (write segments id)
matf = [dr0,'mat/seg_info.mat']; (write semnets info)

& shpf = [dr0,'vector/startpoint.shp']; (read seedpoint shapefile if you have)
% shpf_all = [dr0,'vector/str_dissolve.shp']; (read dissolve shapefile if you have)
pxsz = fdrinfo.pxszx; (read pixel size)
```

  * Stream Network.<br>
The synthetic stream network can be adjusted by customizing the size of the minimum catchment:<br>

```Matlab
mi2px = 5280*.3048/pxsz; % horizontal units assumed to be meters (otherwise modify as needed)
seglen = round(5*mi2px); % 5 miles maximum segment length
sqmi2px = (5280^2)*(.3048^2)/(pxsz^2); % horizontal units assumed to be meters (otherwise modify as needed)
facthr = round(25*sqmi2px); % 25 sq. mi minimum catchment size (min for segment break)
strthr = round(70*sqmi2px); % 70 sq. mi minimum catchment size (min for defining a stream)
```

  * Seed Points.<br>
The FLDPLN stream network can also be specified using seed points:<br>

```Matlab
create_segdb_facthr_maxlen(strthr,facthr,seglen,fdrf,facf,segf,matf) % With no seedpts
% create_segdb_facthr_maxlen_seedpts(shpf,facthr,seglen,fdrf,facf,segf,matf) % With seedpts
```

  * Parallelization.<br>
Parallel execution currently uses the parfeval method from the Matlab parallel package. If you want to run using a serial for loop (required for execution in Octave), comment the following code:<br>

```Matlab
*  Worker Pool
======================================================================================================
if(1>0)
  g = gcp('nocreate');
  if(~isempty(g))
    p = gcp;
  else
    p  = parpool;
  end
  for j = 1:num
    pf(j) = parfeval(p,@fldpln_model_v5ram,1,seg_list(j),inp,fildat,fdrdat,fldmn,fldmx,dh,mxht,1,bg);
  end
% wait(pf)
  pfresults = cell(num,1);
  for idx = 1:num
% fetchNext blocks until the next results are available.
    [completedIdx,value] = fetchNext(pf);
    pfresults{completedIdx} = value;
    fprintf('Segment %d completed/n', completedIdx);
  end
% delete(gcp('nocreate'))
end
return;
======================================================================================================
```

and uncomment the following lines
```Matlab
=================================================================================
% parfor j = 1:num
for j = 1:num
  fldpln_model_v5ram(seg_list(j),inp,fildat,fdrdat,fldmn,fldmx,dh,mxht,1,bg);
end
return;
=================================================================================
```

<br>
<br>
<br>

**[3] create_fldpln_map_cks_Tuscaloosademo.m<br>**

  * Initialize file names.<br>

```Matlab
        fdrf = [dr0,'bil/Fdr.bil']; (read flow direction)
        segdr = [dr0,'segment_files/']; (read segment data in segment_folder)
        shpf = [dr0,'vector/str_segid_polyline_Dissolve.shp']; (read Dissolve shapefile created in ArcGIS)
 ```   

  * Read spatial parameters from FDR.<br>

```Matlab
outf = [dr0,'bil/fldpln_map_h',hs,'_dh',dhs,'_',aoi,'.bil']; (FLDPLN output)
fldmn = 0.01; (minimum depth)

seg_list = sort(readshapefield(shpf,{'grid_code'})); (customize grid_code if needed)

custom_floodplain_map(r,c,h,dh,fldmn,seg_list,h,segdr,outf,fdrf)
return;
hts = h*ones(num,1);
custom_floodplain_map(r,c,h,dh,fldmn,seg_list,hts,segdr,outf,fdrf)
return;
```

<br>
<br>
<br>

**[4] Depth to Flood (DTF) map in ArcGIS<br>**
  
  * After completion of [3] create_fldpln_map_cks_Tuscaloosademo.m, you will have a `fldpln_map_h3_dh0p5_Tuscaloosa.bil` in the `bil` folder, and it will show you a DTF map in ArcGIS.<br>
  
![Result](https://github.com/jameshalgren/fldpln/assets/122640336/f71e7aed-0c35-4fe1-86be-8d5683df90f7)

<br>
<br>
<br>
