#' Draw an area-proportional Venn diagram of 2 or 3 circles
#'
#' This function creates an area-proportional Venn diagram of 2 or 3 circles, based on lists of (biological) identifiers.
#' It requires three parameters: input lists X, Y and Z. For a 2-circle Venn diagram, one of these lists
#' should be left empty. Duplicate identifiers are removed automatically, and a mapping from Entrez and/or
#' Affymetrix to Ensembl IDs is available. BioVenn is case-sensitive. In SVG mode, text and numbers can be dragged and dropped.
#'
#' When using a BioVenn diagram for a publication, please cite:
#' BioVenn - a web application for the comparison and visualization of biological lists using area-proportional Venn diagrams
#' T. Hulsen, J. de Vlieg and W. Alkema, BMC Genomics 2008, 9 (1): 488
#'
#' @param list_x (Required) List with IDs from dataset X
#' @param list_y (Required) List with IDs from dataset Y
#' @param list_z (Required) List with IDs from dataset Z
#' @param title (Optional) The title of the Venn diagram (default is "BioVenn")
#' @param t_f (Optional) The font of the main title (default is "serif")
#' @param t_fb (Optional) The font "face" of the main title (1=plain, 2=bold, 3=italic, 4=bold-italic; default is 2)
#' @param t_s (Optional) The size of the main title (cex; relative to the standard size; default is 1.5)
#' @param t_c (Optional) The colour of the main title (default is "black")
#' @param subtitle (Optional) The subtitle of the Venn diagram (default is "(C) 2007-2020 Tim Hulsen")
#' @param st_f (Optional) The font of the subtitle (default is "serif")
#' @param st_fb (Optional) The font "face" of the subtitle (1=plain, 2=bold, 3=italic, 4=bold-italic; default is 2)
#' @param st_s (Optional) The size of the subtitle (cex; relative to the standard size; default is 1.2)
#' @param st_c (Optional) The colour of the subtitle (default is "black")
#' @param xtitle (Optional) The X title of the Venn diagram (default is "ID set X")
#' @param xt_f (Optional) The font of the X title (default is "serif")
#' @param xt_fb (Optional) The font "face" of the X title (1=plain, 2=bold, 3=italic, 4=bold-italic; default is 2)
#' @param xt_s (Optional) The size of the X title (cex; relative to the standard size; default is 1)
#' @param xt_c (Optional) The colour of the X title (default is "black")
#' @param ytitle (Optional) The Y title of the Venn diagram (default is "ID set Y")
#' @param yt_f (Optional) The font of the Y title (default is "serif")
#' @param yt_fb (Optional) The font "face" of the Y title (1=plain, 2=bold, 3=italic, 4=bold-italic; default is 2)
#' @param yt_s (Optional) The size of the Y title (cex; relative to the standard size; default is 1)
#' @param yt_c (Optional) The colour of the Y title (default is "black")
#' @param ztitle (Optional) The Z title of the Venn diagram (default is "ID set Z")
#' @param zt_f (Optional) The font of the Z title (default is "serif")
#' @param zt_fb (Optional) The font "face" of the Z title (1=plain, 2=bold, 3=italic, 4=bold-italic; default is 2)
#' @param zt_s (Optional) The size of the Z title (cex; relative to the standard size; default is 1)
#' @param zt_c (Optional) The colour of the Z title (default is "black")
#' @param nrtype (Optional) The type of the numbers to be displayed: absolute (abs) numbers or percentages (pct) (default is "abs")
#' @param nr_f (Optional) The font of the numbers (default is "serif")
#' @param nr_fb (Optional) The font "face" of the numbers (1=plain, 2=bold, 3=italic, 4=bold-italic; default is 2)
#' @param nr_s (Optional) The size of the numbers (cex; relative to the standard size; default is 1)
#' @param nr_c (Optional) The colour of the numbers (default is "black")
#' @param x_c (Optional) The colour of the X circle (default is "red")
#' @param y_c (Optional) The colour of the X circle (default is "green")
#' @param z_c (Optional) The colour of the X circle (default is "blue")
#' @param bg_c (Optional) The background colour (default is "white")
#' @param width (Optional) The width of the output file (in pixels for BMP/JPEG/PNG/TIF or in centiinch for PDF/SVG; default is 1000)
#' @param height (Optional) The height of the output file (in pixels for BMP/JPEG/PNG/TIF or in centiinch for PDF/SVG; default is 1000)
#' @param output (Optional) Output format: "bmp","jpg","pdf","png","svg" or "tif" (anything else writes to the screen; default is "screen")
#' @param filename (Optional) The name of the output file (default is "biovenn" + extension of the selected output format)
#' @param map2ens (Optional) Map from Entrez or Affymetrix IDs to Ensembl IDs (default is FALSE)
#' @return An image of the Venn diagram is generated in the desired output format.
#' @return Also returns an object with thirteen lists: X, Y, Z, X only, Y only, Z only, XY, XZ, YZ, XY only, XZ only, YZ only, XYZ.
#' @import biomaRt graphics grDevices plotrix svglite
#' @examples
#' list_x <- c("1007_s_at","1053_at","117_at","121_at","1255_g_at","1294_at")
#' list_y <- c("1255_g_at","1294_at","1316_at","1320_at","1405_i_at")
#' list_z <- c("1007_s_at","1405_i_at","1255_g_at","1431_at","1438_at","1487_at","1494_f_at")
#' biovenn <- draw.venn(list_x, list_y, list_z, subtitle="Example diagram", nrtype="abs")
#' @export
draw.venn <- function(list_x, list_y, list_z, title="BioVenn", t_f="serif", t_fb=2, t_s=1.5, t_c="black", subtitle="(C) 2007-2020 Tim Hulsen", st_f="serif", st_fb=2, st_s=1.2, st_c="black", xtitle="ID Set X", xt_f="serif", xt_fb=2, xt_s=1, xt_c="black", ytitle="ID Set Y", yt_f="serif", yt_fb=2, yt_s=1, yt_c="black", ztitle="ID Set Z", zt_f="serif", zt_fb=2, zt_s=1, zt_c="black", nrtype="abs", nr_f="serif", nr_fb=2, nr_s=1, nr_c="black", x_c="red", y_c="green", z_c="blue", bg_c="white", width=1000, height=1000, output="screen", filename=NULL, map2ens=FALSE){
  # Make input lists unique
  list_x <- unique(list_x)
  list_y <- unique(list_y)
  list_z <- unique(list_z)

  # Convert to Ensembl IDs
  if(map2ens)
  {
    mart <- biomaRt::useMart(dataset="hsapiens_gene_ensembl",biomart="ensembl")
    if(length(list_x)>0)
    {
      list_x_1 <- biomaRt::select(mart, keys=list_x, columns=c("ensembl_gene_id"),keytype="affy_hg_u133a")$ensembl_gene_id
      list_x_2 <- biomaRt::select(mart, keys=list_x, columns=c("ensembl_gene_id"),keytype="entrezgene_id")$ensembl_gene_id
      list_x <- unique(c(list_x_1,list_x_2))
    }
    if(length(list_y)>0)
    {
      list_y_1 <- biomaRt::select(mart, keys=list_y, columns=c("ensembl_gene_id"),keytype="affy_hg_u133a")$ensembl_gene_id
      list_y_2 <- biomaRt::select(mart, keys=list_y, columns=c("ensembl_gene_id"),keytype="entrezgene_id")$ensembl_gene_id
      list_y <- unique(c(list_y_1,list_y_2))
    }
    if(length(list_z)>0)
    {
      list_z_1 <- biomaRt::select(mart, keys=list_z, columns=c("ensembl_gene_id"),keytype="affy_hg_u133a")$ensembl_gene_id
      list_z_2 <- biomaRt::select(mart, keys=list_z, columns=c("ensembl_gene_id"),keytype="entrezgene_id")$ensembl_gene_id
      list_z <- unique(c(list_z_1,list_z_2))
    }
  }

  # Generate lists and calculate numbers
  x <- length(list_x)
  y <- length(list_y)
  z <- length(list_z)
  list_xy <- intersect(list_x, list_y)
  xy <- length(list_xy)
  list_xz <- intersect(list_x, list_z)
  xz <- length(list_xz)
  list_yz <- intersect(list_y, list_z)
  yz <- length(list_yz)
  list_xyz <- intersect(list_xy, list_z)
  xyz <- length(list_xyz)
  list_xy_only <- setdiff(list_xy, list_xyz)
  xy_only <- length(list_xy_only)
  list_xz_only <- setdiff(list_xz, list_xyz)
  xz_only <- length(list_xz_only)
  list_yz_only <- setdiff(list_yz, list_xyz)
  yz_only <- length(list_yz_only)
  list_x_only <- setdiff(list_x, c(list_xy, list_xz))
  x_only <- length(list_x_only)
  list_y_only <- setdiff(list_y, c(list_xy,list_yz))
  y_only <- length(list_y_only)
  list_z_only <- setdiff(list_z, c(list_xz,list_yz))
  z_only <- length(list_z_only)

  # Print numerical output
  print(paste("x total:",x))
  print(paste("y total:",y))
  print(paste("z total:",z))
  print(paste("x only:",x_only))
  print(paste("y only:",y_only))
  print(paste("z only:",z_only))
  print(paste("x-y total overlap:",xy))
  print(paste("x-z total overlap:",xz))
  print(paste("y-z total overlap:",yz))
  print(paste("x-y only overlap:",xy_only))
  print(paste("x-z only overlap:",xz_only))
  print(paste("y-z only overlap:",yz_only))
  print(paste("x-y-z overlap:",xyz))

  # Define sq function
  sq <- function(nr)
  {
    nr=nr^2
    return(nr)
  }

  # Define sqr function
  sqr <- function(nr)
  {
    nr=sqrt(round(abs(nr)))
    return(nr)
  }

  # Define arccos function
  arccos <- function(nr)
  {
    nr=acos(round(nr,5))
    return(nr)
  }

  # Set width and height of plotting area
  width_p=1000
  height_p=1000

  # Amplification
  amp=100000/(x+y+z-xy-xz-yz+xyz)
  x_text=x
  x=x*amp
  y_text=y
  y=y*amp
  z_text=z
  z=z*amp
  xy_text=xy
  xy=xy*amp
  xz_text=xz
  xz=xz*amp
  yz_text=yz
  yz=yz*amp
  xyz_text=xyz
  xyz=xyz*amp
  total=x+y+z-xy-xz-yz+xyz
  total_text=x_text+y_text+z_text-xy_text-xz_text-yz_text+xyz_text

  # Radius calculation
  x_r=sqr(x/pi)
  y_r=sqr(y/pi)
  z_r=sqr(z/pi)

  # Distance calculation
  xy_d=x_r+y_r
  if(x&&y)
  {
    while(xy>sq(x_r)*arccos((sq(xy_d)+sq(x_r)-sq(y_r))/(2*xy_d*x_r))+sq(y_r)*arccos((sq(xy_d)+sq(y_r)-sq(x_r))/(2*xy_d*y_r))-0.5*sqr(round((xy_d+x_r+y_r)*(xy_d+x_r-y_r)*(xy_d-x_r+y_r)*(-xy_d+x_r+y_r),5)))
    {
      xy_d=xy_d-min(x_r,y_r)/1000.0
    }
  }
  xz_d=x_r+z_r
  if(x&&z)
  {
    while(xz>sq(x_r)*arccos((sq(xz_d)+sq(x_r)-sq(z_r))/(2*xz_d*x_r))+sq(z_r)*arccos((sq(xz_d)+sq(z_r)-sq(x_r))/(2*xz_d*z_r))-0.5*sqr(round((xz_d+x_r+z_r)*(xz_d+x_r-z_r)*(xz_d-x_r+z_r)*(-xz_d+x_r+z_r),5)))
    {
      xz_d=xz_d-min(x_r,z_r)/1000.0
    }
  }
  yz_d=y_r+z_r
  if(y&&z)
  {
    while(yz>sq(y_r)*arccos((sq(yz_d)+sq(y_r)-sq(z_r))/(2*yz_d*y_r))+sq(z_r)*arccos((sq(yz_d)+sq(z_r)-sq(y_r))/(2*yz_d*z_r))-0.5*sqr(round((yz_d+y_r+z_r)*(yz_d+y_r-z_r)*(yz_d-y_r+z_r)*(-yz_d+y_r+z_r),5)))
    {
      yz_d=yz_d-min(y_r,z_r)/1000.0
    }
  }
  if(xy_d>xz_d+yz_d){xy_d=xz_d+yz_d;}
  if(xz_d>xy_d+yz_d){xz_d=xy_d+yz_d;}
  if(yz_d>xy_d+xz_d){yz_d=xy_d+xz_d;}

  # Angle calculation
  x_a=arccos((sq(xy_d)+sq(xz_d)-sq(yz_d))/(2*xy_d*xz_d))
  y_a=arccos((sq(xy_d)+sq(yz_d)-sq(xz_d))/(2*xy_d*yz_d))
  z_a=arccos((sq(xz_d)+sq(yz_d)-sq(xy_d))/(2*xz_d*yz_d))
  x_yz=xz_d*sin(z_a)
  y_yz=xy_d*cos(y_a)

  # PPU calculation
  width_h=max(y_r+y_yz,x_r,z_r-yz_d+y_yz)+max(x_r,y_r-y_yz,z_r+yz_d-y_yz)
  ppu_h=width_p/width_h
  width_v=max(x_r+x_yz,y_r,z_r)+max(y_r,z_r,x_r-x_yz)
  ppu_v=height_p/width_v
  ppu=min(ppu_h,ppu_v)

  # Circle center calculation
  x_h=max(x_r,y_r+y_yz,z_r-yz_d+y_yz)
  x_v=max(x_r,y_r-x_yz,z_r-x_yz)
  y_h=max(x_r-y_yz,y_r,z_r-yz_d)
  y_v=max(x_r+x_yz,y_r,z_r)
  z_h=max(x_r+yz_d-y_yz,y_r+yz_d,z_r)
  z_v=max(x_r+x_yz,y_r,z_r)

  # Calculate intersection points X-Y (first inner, then outer)
  xy_i_h_part1=(x_h+y_h)/2+((y_h-x_h)*(sq(x_r)-sq(y_r)))/(2*sq(xy_d))
  xy_i_v_part1=(x_v+y_v)/2+((y_v-x_v)*(sq(x_r)-sq(y_r)))/(2*sq(xy_d))
  xy_i_h_part2=2*((x_v-y_v)/sq(xy_d))*sqr((xy_d+x_r+y_r)*(xy_d+x_r-y_r)*(xy_d-x_r+y_r)*(-xy_d+x_r+y_r))/4
  xy_i_v_part2=2*((x_h-y_h)/sq(xy_d))*sqr((xy_d+x_r+y_r)*(xy_d+x_r-y_r)*(xy_d-x_r+y_r)*(-xy_d+x_r+y_r))/4
  xy_i1_h=xy_i_h_part1-xy_i_h_part2
  xy_i1_v=xy_i_v_part1+xy_i_v_part2
  xy_i2_h=xy_i_h_part1+xy_i_h_part2
  xy_i2_v=xy_i_v_part1-xy_i_v_part2

  # Calculate intersection points X-Z (first inner, then outer)
  xz_i_h_part1=(x_h+z_h)/2+((z_h-x_h)*(sq(x_r)-sq(z_r)))/(2*sq(xz_d))
  xz_i_v_part1=(x_v+z_v)/2+((z_v-x_v)*(sq(x_r)-sq(z_r)))/(2*sq(xz_d))
  xz_i_h_part2=2*((x_v-z_v)/sq(xz_d))*sqr((xz_d+x_r+z_r)*(xz_d+x_r-z_r)*(xz_d-x_r+z_r)*(-xz_d+x_r+z_r))/4
  xz_i_v_part2=2*((x_h-z_h)/sq(xz_d))*sqr((xz_d+x_r+z_r)*(xz_d+x_r-z_r)*(xz_d-x_r+z_r)*(-xz_d+x_r+z_r))/4
  xz_i1_h=xz_i_h_part1+xz_i_h_part2
  xz_i1_v=xz_i_v_part1-xz_i_v_part2
  xz_i2_h=xz_i_h_part1-xz_i_h_part2
  xz_i2_v=xz_i_v_part1+xz_i_v_part2

  # Calculate intersection points Y-Z (first inner, then outer)
  yz_i_h_part1=(y_h+z_h)/2+((z_h-y_h)*(sq(y_r)-sq(z_r)))/(2*sq(yz_d))
  yz_i_v_part1=(y_v+z_v)/2+((z_v-y_v)*(sq(y_r)-sq(z_r)))/(2*sq(yz_d))
  yz_i_h_part2=2*((y_v-z_v)/sq(yz_d))*sqr((yz_d+y_r+z_r)*(yz_d+y_r-z_r)*(yz_d-y_r+z_r)*(-yz_d+y_r+z_r))/4
  yz_i_v_part2=2*((y_h-z_h)/sq(yz_d))*sqr((yz_d+y_r+z_r)*(yz_d+y_r-z_r)*(yz_d-y_r+z_r)*(-yz_d+y_r+z_r))/4
  yz_i1_h=yz_i_h_part1-yz_i_h_part2
  yz_i1_v=yz_i_v_part1+yz_i_v_part2
  yz_i2_h=yz_i_h_part1+yz_i_h_part2
  yz_i2_v=yz_i_v_part1-yz_i_v_part2

  # Number fill point calculation of overlaps
  xy_f_h=(xy_i2_h+xz_i1_h+yz_i1_h)/3
  xy_f_v=(xy_i2_v+xz_i1_v+yz_i1_v)/3
  xz_f_h=(xy_i1_h+xz_i2_h+yz_i1_h)/3
  xz_f_v=(xy_i1_v+xz_i2_v+yz_i1_v)/3
  yz_f_h=(xy_i1_h+xz_i1_h+yz_i2_h)/3
  yz_f_v=(xy_i1_v+xz_i1_v+yz_i2_v)/3
  xyz_f_h=(xy_i1_h+xz_i1_h+yz_i1_h)/3
  xyz_f_v=(xy_i1_v+xz_i1_v+yz_i1_v)/3

  # Number fill point calculation of X
  # For XYZ diagrams
  if(x&&y&&z)
  {
    xyz_yz_i1=sqr(sq(xyz_f_h-yz_i1_h)+sq(xyz_f_v-yz_i1_v))
    x_ratio_h=(xyz_f_h-yz_i1_h)/xyz_yz_i1
    x_ratio_v=(xyz_f_v-yz_i1_v)/xyz_yz_i1
    x_out_h=x_h-x_r*x_ratio_h
    x_out_v=x_v-x_r*x_ratio_v
    x_f_h=(x_out_h+yz_i1_h)/2
    x_f_v=(x_out_v+yz_i1_v)/2
  }
  # For XY diagrams
  else if(x&&y&&!z)
  {
    xy_f_h=(xy_i1_h+xy_i2_h)/2
    xy_f_v=(xy_i1_v+xy_i2_v)/2
    x_in_h=y_h+cos(y_a)*y_r
    x_in_v=y_v-sin(y_a)*y_r
    x_out_h=x_h+cos(y_a)*x_r
    x_out_v=x_v-sin(y_a)*x_r
    x_f_h=(x_out_h+x_in_h)/2
    x_f_v=(x_out_v+x_in_v)/2
  }
  # For XZ diagrams
  else if(x&&!y&&z)
  {
    xz_f_h=(xz_i1_h+xz_i2_h)/2
    xz_f_v=(xz_i1_v+xz_i2_v)/2
    x_in_h=z_h-cos(z_a)*z_r
    x_in_v=z_v-sin(z_a)*z_r
    x_out_h=x_h-cos(z_a)*x_r
    x_out_v=x_v-sin(z_a)*x_r
    x_f_h=(x_out_h+x_in_h)/2
    x_f_v=(x_out_v+x_in_v)/2
  }

  # Number fill point calculation of Y
  # For XYZ diagrams
  if(x&&y&&z)
  {
    xyz_xz_i1=sqr(sq(xyz_f_h-xz_i1_h)+sq(xyz_f_v-xz_i1_v))
    y_ratio_h=(xyz_f_h-xz_i1_h)/xyz_xz_i1
    y_ratio_v=(xyz_f_v-xz_i1_v)/xyz_xz_i1
    y_out_h=y_h-y_r*y_ratio_h
    y_out_v=y_v-y_r*y_ratio_v
    y_f_h=(y_out_h+xz_i1_h)/2
    y_f_v=(y_out_v+xz_i1_v)/2
  }
  # For XY diagrams
  else if(x&&y&&!z)
  {
    xy_f_h=(xy_i1_h+xy_i2_h)/2
    xy_f_v=(xy_i1_v+xy_i2_v)/2
    y_in_h=x_h-cos(y_a)*x_r
    y_in_v=x_v+sin(y_a)*x_r
    y_out_h=y_h-cos(y_a)*y_r
    y_out_v=y_v+sin(y_a)*y_r
    y_f_h=(y_out_h+y_in_h)/2
    y_f_v=(y_out_v+y_in_v)/2
  }
  # For YZ diagrams
  else if(!x&&y&&z)
  {
    yz_f_h=(yz_i1_h+yz_i2_h)/2
    yz_f_v=(yz_i1_v+yz_i2_v)/2
    y_in_h=z_h-z_r
    y_in_v=z_v
    y_out_h=y_h-y_r
    y_out_v=y_v
    y_f_h=(y_out_h+y_in_h)/2
    y_f_v=(y_out_v+y_in_v)/2
  }

  # Number fill point calculation of Z
  # For XYZ diagrams
  if(x&&y&&z)
  {
    xyz_xy_i1=sqr(sq(xyz_f_h-xy_i1_h)+sq(xyz_f_v-xy_i1_v))
    z_ratio_h=(xyz_f_h-xy_i1_h)/xyz_xy_i1
    z_ratio_v=(xyz_f_v-xy_i1_v)/xyz_xy_i1
    z_out_h=z_h-z_r*z_ratio_h
    z_out_v=z_v-z_r*z_ratio_v
    z_f_h=(z_out_h+xy_i1_h)/2
    z_f_v=(z_out_v+xy_i1_v)/2
  }
  # For XZ diagrams
  else if(x&&!y&&z)
  {
    xz_f_h=(xz_i1_h+xz_i2_h)/2
    xz_f_v=(xz_i1_v+xz_i2_v)/2
    z_in_h=x_h+cos(z_a)*x_r
    z_in_v=x_v+sin(z_a)*x_r
    z_out_h=z_h+cos(z_a)*z_r
    z_out_v=z_v+sin(z_a)*z_r
    z_f_h=(z_out_h+z_in_h)/2
    z_f_v=(z_out_v+z_in_v)/2
  }
  # For YZ diagrams
  else if(!x&&y&&z)
  {
    yz_f_h=(yz_i1_h+yz_i2_h)/2
    yz_f_v=(yz_i1_v+yz_i2_v)/2
    z_in_h=y_h+z_r
    z_in_v=y_v
    z_out_h=z_h+y_r
    z_out_v=z_v
    z_f_h=(z_out_h+z_in_h)/2
    z_f_v=(z_out_v+z_in_v)/2
  }

  # Number fill point calculation for special cases
  # No X only
  if(x&&!x_only)
  {
    # X is subset of Z
    if(!xy_only)
    {
      xz_f_h=((y_h+y_r)+(x_h+x_r))/2
      xz_f_v=x_v
      yz_f_h=((z_h-z_r)+(x_h-x_r))/2
      yz_f_v=y_v
      xyz_f_h=((x_h-x_r)+(y_h+y_r))/2
      xyz_f_v=(x_v+y_v)/2
      y_f_h=((y_h-y_r)+(z_h-z_r))/2
      y_f_v=y_v
      xyz_xy_i2=sqr(sq(xyz_f_h-xy_i2_h)+sq(xyz_f_v-xy_i2_v))
      z_ratio_h=(xyz_f_h-xy_i2_h)/xyz_xy_i2
      z_ratio_v=(xyz_f_v-xy_i2_v)/xyz_xy_i2
      z_out_h=z_h-z_r*z_ratio_h
      z_out_v=z_v-z_r*z_ratio_v
      z_f_h=(z_out_h+xy_i2_h)/2
      z_f_v=(z_out_v+xy_i2_v)/2
    }
    # X is subset of Y
    else if(!xz_only)
    {
      xy_f_h=((y_h-y_r)+(z_h-z_r))/2
      xy_f_v=x_v
      yz_f_h=((x_h+x_r)+(y_h+y_r))/2
      yz_f_v=z_v
      xyz_f_h=((z_h-z_r)+(x_h+x_r))/2
      xyz_f_v=(x_v+z_v)/2
      xyz_xz_i2=sqr(sq(xyz_f_h-xz_i2_h)+sq(xyz_f_v-xz_i2_v))
      y_ratio_h=(xyz_f_h-xz_i2_h)/xyz_xz_i2
      y_ratio_v=(xyz_f_v-xz_i2_v)/xyz_xz_i2
      y_out_h=y_h-y_r*y_ratio_h
      y_out_v=y_v-y_r*y_ratio_v
      y_f_h=(y_out_h+xz_i2_h)/2
      y_f_v=(y_out_v+xz_i2_v)/2
      z_f_h=((y_h+y_r)+(z_h+z_r))/2
      z_f_v=z_v
    }
  }
  # No Y only
  if(y&&!y_only)
  {
    # Y is subset of Z
    if(!xy_only)
    {
      xz_f_h=((y_h+y_r)+(z_h+z_r))/2
      xz_f_v=x_v
      yz_f_h=((x_h-x_r)+(y_h-y_r))/2
      yz_f_v=y_v
      xyz_f_h=((x_h-x_r)+(y_h+y_r))/2
      xyz_f_v=(x_v+y_v)/2
      x_f_h=((z_h+z_r)+(x_h+x_r))/2
      x_f_v=x_v
      xyz_xy_i2=sqr(sq(xyz_f_h-xy_i2_h)+sq(xyz_f_v-xy_i2_v))
      z_ratio_h=(xyz_f_h-xy_i2_h)/xyz_xy_i2
      z_ratio_v=(xyz_f_v-xy_i2_v)/xyz_xy_i2
      z_out_h=z_h-z_r*z_ratio_h
      z_out_v=z_v-z_r*z_ratio_v
      z_f_h=(z_out_h+xy_i2_h)/2
      z_f_v=(z_out_v+xy_i2_v)/2
    }
    # Y is subset of X
    else if(!yz_only)
    {
      xy_f_h=((y_h-y_r)+(z_h-z_r))/2
      xy_f_v=y_v
      xz_f_h=((y_h+y_r)+(x_h+x_r))/2
      xz_f_v=z_v
      xyz_f_h=((z_h-z_r)+(y_h+y_r))/2
      xyz_f_v=(y_v+z_v)/2
      xyz_yz_i1=sqr(sq(xyz_f_h-yz_i1_h)+sq(xyz_f_v-yz_i1_v))
      x_ratio_h=(xyz_f_h-yz_i1_h)/xyz_yz_i1
      x_ratio_v=(xyz_f_v-yz_i1_v)/xyz_yz_i1
      x_out_h=x_h-x_r*x_ratio_h
      x_out_v=x_v-x_r*x_ratio_v
      x_f_h=(x_out_h+yz_i1_h)/2
      x_f_v=(x_out_v+yz_i1_v)/2
      z_f_h=((x_h+x_r)+(z_h+z_r))/2
      z_f_v=y_v
    }
  }
  # No Z only
  if(z&&!z_only)
  {
    # Z is subset of Y
    if(!xz_only)
    {
      xy_f_h=((y_h-y_r)+(z_h-z_r))/2
      xy_f_v=x_v
      yz_f_h=((x_h+x_r)+(z_h+z_r))/2
      yz_f_v=z_v
      xyz_f_h=((z_h-z_r)+(x_h+x_r))/2
      xyz_f_v=(x_v+z_v)/2
      xyz_xz_i2=sqr(sq(xyz_f_h-xz_i2_h)+sq(xyz_f_v-xz_i2_v))
      y_ratio_h=(xyz_f_h-xz_i2_h)/xyz_xz_i2
      y_ratio_v=(xyz_f_v-xz_i2_v)/xyz_xz_i2
      y_out_h=y_h-y_r*y_ratio_h
      y_out_v=y_v-y_r*y_ratio_v
      y_f_h=(y_out_h+xz_i2_h)/2
      y_f_v=(y_out_v+xz_i2_v)/2
      x_f_h=((x_h-x_r)+(y_h-y_r))/2
      x_f_v=x_v
    }
    # Z is subset of X
    else if(!yz_only)
    {
      xy_f_h=((z_h-z_r)+(x_h-x_r))/2
      xy_f_v=y_v
      xz_f_h=((y_h+y_r)+(x_h+x_r))/2
      xz_f_v=x_v
      xyz_f_h=((z_h-z_r)+(y_h+y_r))/2
      xyz_f_v=(y_v+z_v)/2
      xyz_yz_i1=sqr(sq(xyz_f_h-yz_i1_h)+sq(xyz_f_v-yz_i1_v))
      x_ratio_h=(xyz_f_h-yz_i1_h)/xyz_yz_i1
      x_ratio_v=(xyz_f_v-yz_i1_v)/xyz_yz_i1
      x_out_h=x_h-x_r*x_ratio_h
      x_out_v=x_v-x_r*x_ratio_v
      x_f_h=(x_out_h+yz_i1_h)/2
      x_f_v=(x_out_v+yz_i1_v)/2
      y_f_h=((y_h-y_r)+(x_h-x_r))/2
      y_f_v=y_v
    }
  }

  # Output to file or screen
  if(output=="bmp")
  {
    if(is.null(filename))
    {
      filename="biovenn.bmp"
    }
    grDevices::bmp(filename,width=width,height=height,units="px")
  }
  else if(output=="jpg")
  {
    if(is.null(filename))
    {
      filename="biovenn.jpg"
    }
    grDevices::jpeg(filename,width=width,height=height,units="px")
  }
  else if(output=="pdf")
  {
    if(is.null(filename))
    {
      filename="biovenn.pdf"
    }
    grDevices::pdf(filename,width=width/100,height=height/100)
  }
  else if(output=="png")
  {
    if(is.null(filename))
    {
      filename="biovenn.png"
    }
    grDevices::png(filename,width=width,height=height,units="px")
  }
  else if(output=="svg")
  {
    if(is.null(filename))
    {
      filename="biovenn.svg"
    }
    svglite::svglite("biovenn_temp.svg",width=width/100,height=height/100)
  }
  else if(output=="tif")
  {
    if(is.null(filename))
    {
      filename="biovenn.tif"
    }
    grDevices::tiff(filename,width=width,height=height,units="px")
  }

  # Draw circles
  opar<-graphics::par(no.readonly=TRUE)
  on.exit(graphics::par(opar))
  graphics::par(pty="s",bg=bg_c)
  graphics::plot(0,type="n",axes=FALSE,xlim=c(0,width_p),ylim=c(height_p,0),xlab="",ylab="",xaxt="none",yaxt="none")
  graphics::par(family=t_f)
  graphics::title(main=title,line=1,font.main=t_fb,cex.main=t_s,col.main=t_c)
  graphics::par(family=st_f)
  graphics::title(sub=subtitle,line=1,font.sub=st_fb,cex.sub=st_s,col.sub=st_c)
  plotrix::draw.circle(ppu*x_h,ppu*x_v,ppu*x_r,lty=0,col=grDevices::rgb(grDevices::col2rgb(x_c)[,1][1],grDevices::col2rgb(x_c)[,1][2],grDevices::col2rgb(x_c)[,1][3],maxColorValue=255,alpha=128))
  plotrix::draw.circle(ppu*y_h,ppu*y_v,ppu*y_r,lty=0,col=grDevices::rgb(grDevices::col2rgb(y_c)[,1][1],grDevices::col2rgb(y_c)[,1][2],grDevices::col2rgb(y_c)[,1][3],maxColorValue=255,alpha=128))
  plotrix::draw.circle(ppu*z_h,ppu*z_v,ppu*z_r,lty=0,col=grDevices::rgb(grDevices::col2rgb(z_c)[,1][1],grDevices::col2rgb(z_c)[,1][2],grDevices::col2rgb(z_c)[,1][3],maxColorValue=255,alpha=128))

  # Print numbers
  if(length(nrtype)>0)
  {
    if(nrtype=="abs")
    {
      if(x_only)
      {
        graphics::text(ppu*x_f_h-nr_s*0.3*length(x_only),ppu*x_f_v,x_only,col=nr_c,family=nr_f,font=nr_fb,cex=nr_s)
      }
      if(y_only)
      {
        graphics::text(ppu*y_f_h-nr_s*0.3*length(y_only),ppu*y_f_v,y_only,col=nr_c,family=nr_f,font=nr_fb,cex=nr_s)
      }
      if(z_only)
      {
        graphics::text(ppu*z_f_h-nr_s*0.3*length(z_only),ppu*z_f_v,z_only,col=nr_c,family=nr_f,font=nr_fb,cex=nr_s)
      }
      if(xy_only)
      {
        graphics::text(ppu*xy_f_h-nr_s*0.3*length(xy_only),ppu*xy_f_v,xy_only,col=nr_c,family=nr_f,font=nr_fb,cex=nr_s)
      }
      if(xz_only)
      {
        graphics::text(ppu*xz_f_h-nr_s*0.3*length(xz_only),ppu*xz_f_v,xz_only,col=nr_c,family=nr_f,font=nr_fb,cex=nr_s)
      }
      if(yz_only)
      {
        graphics::text(ppu*yz_f_h-nr_s*0.3*length(yz_only),ppu*yz_f_v,yz_only,col=nr_c,family=nr_f,font=nr_fb,cex=nr_s)
      }
      if(xyz)
      {
        graphics::text(ppu*xyz_f_h-nr_s*0.3*length(xyz_text),ppu*xyz_f_v,xyz_text,col=nr_c,family=nr_f,font=nr_fb,cex=nr_s)
      }
    }
    else if(nrtype=="pct")
    {
      if(x_only)
      {
        graphics::text(ppu*x_f_h-nr_s*0.3*length(paste0(round(x_only/total_text*100,2),"%")),ppu*x_f_v,paste0(round(x_only/total_text*100,2),"%"),col=nr_c,family=nr_f,font=nr_fb,cex=nr_s)
      }
      if(y_only)
      {
        graphics::text(ppu*y_f_h-nr_s*0.3*length(paste0(round(y_only/total_text*100,2),"%")),ppu*y_f_v,paste0(round(y_only/total_text*100,2),"%"),col=nr_c,family=nr_f,font=nr_fb,cex=nr_s)
      }
      if(z_only)
      {
        graphics::text(ppu*z_f_h-nr_s*0.3*length(paste0(round(z_only/total_text*100,2),"%")),ppu*z_f_v,paste0(round(z_only/total_text*100,2),"%"),col=nr_c,family=nr_f,font=nr_fb,cex=nr_s)
      }
      if(xy_only)
      {
        graphics::text(ppu*xy_f_h-nr_s*0.3*length(paste0(round(xy_only/total_text*100,2),"%")),ppu*xy_f_v,paste0(round(xy_only/total_text*100,2),"%"),col=nr_c,family=nr_f,font=nr_fb,cex=nr_s)
      }
      if(xz_only)
      {
        graphics::text(ppu*xz_f_h-nr_s*0.3*length(paste0(round(xz_only/total_text*100,2),"%")),ppu*xz_f_v,paste0(round(xz_only/total_text*100,2),"%"),col=nr_c,family=nr_f,font=nr_fb,cex=nr_s)
      }
      if(yz_only)
      {
        graphics::text(ppu*yz_f_h-nr_s*0.3*length(paste0(round(yz_only/total_text*100,2),"%")),ppu*yz_f_v,paste0(round(yz_only/total_text*100,2),"%"),col=nr_c,family=nr_f,font=nr_fb,cex=nr_s)
      }
      if(xyz)
      {
        graphics::text(ppu*xyz_f_h-nr_s*0.3*length(paste0(round(xyz_text/total_text*100,2),"%")),ppu*xyz_f_v,paste0(round(xyz_text/total_text*100,2),"%"),col=nr_c,family=nr_f,font=nr_fb,cex=nr_s)
      }
    }
  }

  # Print texts
  if(x)
  {
    graphics::text(ppu*x_h-xt_s*0.3*length(xtitle),ppu*x_v,xtitle,col=xt_c,family=xt_f,font=xt_fb,cex=xt_s)
  }
  if(y)
  {
    graphics::text(ppu*y_h-yt_s*0.3*length(ytitle),ppu*y_v,ytitle,col=yt_c,family=yt_f,font=yt_fb,cex=yt_s)
  }
  if(z)
  {
    graphics::text(ppu*z_h-zt_s*0.3*length(ztitle),ppu*z_v,ztitle,col=zt_c,family=zt_f,font=zt_fb,cex=zt_s)
  }

  # Write to file
  if(output %in% c("bmp","jpg","pdf","png","svg","tif"))
  {
    grDevices::dev.off()
  }

  # Create drag-and-drop functionality for SVG file
  if(output=="svg")
  {
    svg_temp <- file("biovenn_temp.svg", "r")
    svg <- file(filename, "w")
    id=1
    while (length(oneLine <- readLines(svg_temp, n=1, warn=FALSE)) > 0) {
      if(substr(oneLine,1,4)=="<svg")
      {
        oneLine=sub("viewBox='0 0 (\\d*\\.\\d*) (\\d*\\.\\d*)'","viewBox='0 0 \\1 \\2' height='\\1' width='\\2'",oneLine)
      }
      if(substr(oneLine,1,5)=="<rect")
      {
        oneLine=sub("<rect","<script>
<![CDATA[
var Root=document.documentElement
standardize(Root)
function standardize(R){
  var Attr={
    'onmouseup':'add(evt)',
    'onmousedown':'grab(evt)',
    'onmousemove':null
  }
  assignAttr(R,Attr)
}
function grab(evt){
  var O=evt.target
  var Attr={
    'onmousemove':'slide(evt,\"'+O.id+'\")',
    'onmouseup':'standardize(Root)'
  }
  assignAttr(Root,Attr)
}
function slide(evt,id){
  if(id!='rect'&&id!='polygon'){
    var o=document.getElementById(id)
    o.setAttributeNS(null, 'x', evt.clientX)
    o.setAttributeNS(null, 'y', evt.clientY)
  }
}
function assignAttr(O,A){
  for (i in A) O.setAttributeNS(null,i, A[i])
}
]]>
</script>
<rect id='rect'",oneLine)
      }
      else if(substr(oneLine,1,5)=="<text")
      {
        oneLine=paste0(substr(oneLine,1,5)," id='t",id,"'",substr(oneLine,6,nchar(oneLine)))
        oneLine=sub("style='","style='cursor:move;",oneLine)
        id=id+1
      }
      else if(substr(oneLine,61,65)=="<text")
      {
        oneLine=paste0(substr(oneLine,1,65)," id='t",id,"'",substr(oneLine,66,nchar(oneLine)))
        oneLine=sub("style='","style='cursor:move;",oneLine)
        id=id+1
      }
      else if(substr(oneLine,1,8)=="<polygon")
      {
        oneLine=paste0(substr(oneLine,1,8)," id='polygon'",substr(oneLine,9,nchar(oneLine)))
      }
      write(oneLine,svg)
    }
    close(svg_temp)
    file.remove("biovenn_temp.svg")
    close(svg)
  }

  # Return lists
  return(list("x"=list_x,"y"=list_y,"z"=list_z,"x_only"=list_x_only,"y_only"=list_y_only,"z_only"=list_z_only,"xy"=list_xy,"xz"=list_xz,"yz"=list_yz,"xy_only"=list_xy_only,"xz_only"=list_xz_only,"yz_only"=list_yz_only,"xyz"=list_xyz))
}
