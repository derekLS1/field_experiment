#README - REQUIREMENTS

#FIELD IMAGE MUST BE A PNG FILE. 
#PLANTS ARE BLACK CIRCLES OF DIAMETER 25 PIXELS.
#THE SOUTHERN END OF THE RULER/WHITEBOARD IS A RED CIRCLE OF DIAMETER 25 PIXELS.
#THE NORTHERN END OF THE RULER/WHITEBOARD IS A RED CIRCLE OF DIAMETER 50 PIXELS. 
#IT IS ASSUMED THE 50 PIXEL DIAMTER RED CIRCLE CORRESPONDS EXACTLY WITH THE GPS COORDINATES.

#THE FILENAME OF THE OVERHEAD PICTURE MUST HAVE WAYPOINT, NUMBER OF EACH PLANT LEFT TO RIGHT, ALL SEPARATED BY "_"
#EXAMPLE --> overhead_picture_filename="wp345_1_2_3_4_5_6.png"

#THE DIRECTORY OF THE OVERHEAD PICTURES MUST HAVE A SLASH ON THE RIGHT SIDE
#EXAMPLE --> "/Users/dlundberg/Documents/abt6/field_experiment/fall_2015/"

#TO ANALYZE THE IMAGE, SOURCE THIS R SCRIPT FILE TO LOAD ALL THE FUNCTIONS. THEN RUN THE FUNCTION "process_field_image"
#THIS FUNCTION REQUIRES AS INPUT 1) overhead_picture_filename and 2) directory
#IT WRITES THE SUMMARY TABLE TO THE SAME DIRECTORY.

#EXAMPLE USAGE:
#process_field_image(overhead_picture_filename="wp345_1_2_3_4_5_6.png", directory="/Users/dlundberg/Documents/abt6/field_experiment/fall_2015/", ruler=306)




##################################################################################################
##################################################################################################
##################################################################################################
##################################################################################################

library(SDMTools)
library(png)

#DEFINE A COUPLE INTERNAL FUNCTIONS
# Calculates the geodesic distance between two points specified by radian latitude/longitude using the
# FUNCTION::: Haversine formula (hf)
gcd.hf <- function(long1, lat1, long2, lat2) {
  R <- 6371 # Earth mean radius [km]
  delta.long <- (long2 - long1)
  delta.lat <- (lat2 - lat1)
  a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
  c <- 2 * asin(min(1,sqrt(a)))
  d = (R * c)*1000
  return(d) # Distance in m
}
#FUNCTION:::  Convert degrees to radians
deg2rad <- function(deg){return(deg*pi/180)}

#FUNCTION::: RECOGNIZE EACH CLUSTER OF POINTS (EACH PLANT OR RULER MARKER)
dotFinder<-function(binarymatrix){
  ConnCompLabel(binarymatrix)->spots
  nrow(spots)->matrixheight
  unique(as.numeric(spots))->spotnames
  spotnames[which(spotnames!=0)]->spotnames
  spotmatrix=matrix(ncol=7, nrow=length(spotnames), data=0)
  #FOR EACH PLANT DOT, MAKES A ROW IN A DATA MATRIX.
  #Row: 1=plant ID, 2=x value of center, 3=y value of center, 4=width, 5=height, 6=nothing, 7=pixels
  for(s in spotnames){
    which(spots==s)->points_in_s
    ceiling(points_in_s/matrixheight)->s_x
    matrixheight-(points_in_s-((s_x-1)*matrixheight))->s_y
    1+max(s_x)-min(s_x)->xlength
    1+max(s_y)-min(s_y)->ylength
    xlength*ylength->size_of_hypothetical_s_square
    if(length(points_in_s)>(size_of_hypothetical_s_square/2)){
      spotmatrix[which(s==spotnames),]<-c(s, round(mean(s_x), 0), round(mean(s_y),0), xlength, ylength, 1,length(points_in_s))
    } 
  }
  return(spotmatrix)
}

#MAIN FUNCTION

process_field_image<-function(overhead_picture_filename="wp345_1_2_3_4_5_6.png",
                                    directory="/Users/dlundberg/Documents/abt6/field_experiment/fall_2015/",
                                    ruler=306){

#LOAD OVERHEAD PICTURE
field<-readPNG(paste(directory, overhead_picture_filename, collapse="", sep=""),native = FALSE)

#MAKE SURE PLENTY OF DIGITS DISPLAYED
options(digits=15)

#PARSES overhead_picture_filename INTO THE PLANT IDS AND WAYPOINT NAME
gsub(".png", "", overhead_picture_filename)->overhead_picture_filename_parts1
gsub(" ", "", overhead_picture_filename_parts1)->overhead_picture_filename_parts1
do.call(rbind, strsplit(overhead_picture_filename_parts1, split="_"))->overhead_picture_filename_parts
overhead_picture_filename_parts[grep("wp", overhead_picture_filename_parts)]->waypoint
as.numeric(gsub("wp", "", waypoint))->waypoint
overhead_picture_filename_parts[grep("wp", overhead_picture_filename_parts, invert=TRUE)]->plantIDs

#CONVERTS THE RGB IMAGE INTO A SINGLE GREYSCALE IMAGE
field[,,1]+field[,,2]+field[,,3]->field
field/3->field

#THRESHOLDS MATRIX OF IMAGE TO SEPARATE RED FROM BLACK
as.matrix(field)->field
1->field[which(field>=0.9)]
0->field[which(field<=0.1)]
0.5->field[which(field>0.2 & field<0.8)]

#SWAPS AROUND COLORS SO 0=WHITE, 1=BLACK, and 0.5=RED
2->field[which(field==1)]
1->field[which(field==0)]
0->field[which(field==2)]

#CREATES TWO 1/0 BINARY MATRIXES FROM THE FILE, ONE FOR THE RULER ("red"), ONE FOR THE PLANTS ("field")
#RULER REMOVED FROM PLANTS PICTURE AND VICE VERSA
red<-field
#remove black from ruler pic
0->red[which(red==1)]
#convert red (0.5) to black (1) in ruler pic
1->red[which(red==0.5)]
#remove ruler (0.5) from plants pic
0->field[which(field==0.5)]

#PLANTS PIC RENAMED AS "BLACK"
field->black

#RECOGNIZE DOTS IN IMAGES AND RECORD THE LOCATION AND DOT SIZE IN A SIMPLE MATRIX
dotFinder(black)->blackmatrix
dotFinder(red)->redmatrix

#FILTER OUT POTENTIAL FALSE DOTS THAT ARE SMALLER THAN 100 TOTAL PIXELS
redmatrix[redmatrix[,7]>100,]->redmatrix
blackmatrix[blackmatrix[,7]>100,]->blackmatrix

#SORT THE REDMATRIX SUCH THAT THE BIG DOT CORRESPONDING TO NORTH IS THE FIRST ROW
redmatrix[order(redmatrix[,7], decreasing=TRUE),]->redmatrix

#SORT THE PLANT DOTS BY THE ORDER IN WHICH THEY APPEAR LEFT TO RIGHT. NAME THEM WITH THE PLANT IDS THAT
#CAME FROM THE FILENAME.
blackmatrix[order(blackmatrix[,2]),]->blackmatrix
if(length(plantIDs)==length(blackmatrix[,1])){
  as.numeric(plantIDs)->blackmatrix[,1]  
}
if(length(plantIDs)!=length(blackmatrix[,1])){
  print("DOT NUMBER DOES NOT MATCH NUMBER OF PLANTS IN FILENAME") 
}

#USE THE RED DOTS OF THE WHITEBOARD TO CALCULATE THE DELTA X AND DELTA Y
(redmatrix[1,2]-redmatrix[2,2])->xruler
(redmatrix[1,3]-redmatrix[2,3])->yruler

#FIND THE ROTATIONAL CORRECTION FOR THE RULER
#if(xruler<1 & yruler<1){}
-atan(xruler/yruler)->radians_rotation
#((pi/2)+radians_rotation)->radians_rotation
round(radians_rotation*180/pi, 2)->degree_rotation
if(yruler<1){
    radians_rotation+pi->radians_rotation
    degree_rotation+180->degree_rotation
}


#COMBINE REDMATRIX AND BLACKMATRIX INTO A SINGLE SUMMARY MATRIX
rbind(redmatrix, blackmatrix)->allmatrix

#APPLY TRIGONOMETRIC ROTATION TO ALL X AND Y VALUES IN REDMATRIX
allmatrix[,2]->x
allmatrix[,3]->y
allmatrix->allmatrix_rotated
cos(radians_rotation)*(x) + sin(radians_rotation)*(y) -> allmatrix_rotated[,2]
-sin(radians_rotation)*(x) + cos(radians_rotation)*(y) -> allmatrix_rotated[,3]

#RESCALE VALUES TO MILLIMETERS BASED ON LENGTH OF RULER AND CENTER PLOT ON "NORTH"
rulerlength=round(sqrt( xruler^2 + yruler^2 ), 4) #finds length of original ruler
ruler/rulerlength->sizing_factor 
allmatrix_rotated->allmatrix_rotated_resized
allmatrix_rotated_resized[,2:7]*sizing_factor->allmatrix_rotated_resized[,2:7]
allmatrix_rotated_resized[1,]->GPS_reference
allmatrix_rotated_resized[,2]-GPS_reference[2]->allmatrix_rotated_resized[,2]
allmatrix_rotated_resized[,3]-GPS_reference[3]->allmatrix_rotated_resized[,3]


gsub("wp", "simple", overhead_picture_filename_parts1)->simplename
     
png(paste(directory, simplename, ".png", sep="", collapse=""), width=600, height=1400)

#GRAPH RESULTS IN THREE PANELS
par(mfrow=c(3,1))
#FIRST PLOT OF UNTRANSFORMED SIMPLIFIED IMAGE
allmatrix->P
plot(P[3:nrow(P),2], P[3:nrow(P),3], col="white", xlim=c(min(P[,2]), max(P[,2])), ylim=c(min(P[,3]), max(P[,3])), xlab="pixels", ylab="pixels", main=paste("waypoint ", waypoint, ", simplified image", sep="", collapse=""))
points(P[1:2,2], P[1:2,3], col="red", type="l")
text(P[1:2,2], P[1:2,3], col="red", labels=c("N", "S"), cex=4)
text(P[3:nrow(P),2], P[3:nrow(P),3], labels=P[3:nrow(P),1])

#SECOND PLOT OF ROTATED IMAGE
allmatrix_rotated->P
plot(P[3:nrow(P),2], P[3:nrow(P),3], col="white", xlim=c(min(P[,2]), max(P[,2])), ylim=c(min(P[,3]), max(P[,3])), xlab="pixels", ylab="pixels", main=paste(degree_rotation, " degree rotation", sep="", collapse="")) 
points(P[1:2,2], P[1:2,3], col="red", type="l")
text(P[1:2,2], P[1:2,3], col="red", labels=c("N", "S"), cex=4)
text(P[3:nrow(P),2], P[3:nrow(P),3], labels=P[3:nrow(P),1])

#THIRD PLOT OF ROTATED AND RESIZED IMAGE
allmatrix_rotated_resized->P
plot(P[3:nrow(P),2], P[3:nrow(P),3], col="white", xlim=c(min(P[,2]), max(P[,2])), ylim=c(min(P[,3]), max(P[,3])), xlab="distance (mm)", ylab="distance (mm)", main=paste(degree_rotation, " degree rotation, rescaled ", round(sizing_factor, 2), sep="", collapse="")) 
points(P[1:2,2], P[1:2,3], col="red", type="l")
text(P[1:2,2], P[1:2,3], col="red", labels=c("N", "S"), cex=4)
text(P[3:nrow(P),2], P[3:nrow(P),3], labels=P[3:nrow(P),1])

dev.off()


#PREPARE SUMMARY MATRIX FOR PRINT
cbind(rep(as.numeric(waypoint), nrow(allmatrix_rotated_resized)), allmatrix_rotated_resized)->allmatrixRR_print
round(allmatrixRR_print, 2)->allmatrixRR_print
allmatrixRR_print[1:2,2]<-c("N","S")
rbind(c("waypoint", "plantID", "horiz_mm_from_GPS", "vert_mm_from_GPS"), allmatrixRR_print[,c(1,2,3,4)])->allmatrixRR_print

#WRITE SUMMARY MATRIX AND PLOT
allmatrix->P
plot(P[3:nrow(P),2], P[3:nrow(P),3], col="white", xlim=c(min(P[,2]), max(P[,2])), ylim=c(min(P[,3]), max(P[,3])), xlab="pixels", ylab="pixels", main=paste("waypoint ", waypoint, ", simplified image", sep="", collapse=""))
points(P[1:2,2], P[1:2,3], col="red", type="l")
text(P[1:2,2], P[1:2,3], col="red", labels=c("N", "S"), cex=4)
text(P[3:nrow(P),2], P[3:nrow(P),3], labels=P[3:nrow(P),1])
write.table(allmatrixRR_print, file = paste(directory, overhead_picture_filename_parts1, ".txt", collapse="", sep=""), sep = "\t",row.names = FALSE, col.names = FALSE, quote=FALSE)

}


#process_field_image(overhead_picture_filename="wp345_1_2_3_4_5_6.png", directory="/Users/dlundberg/Documents/abt6/field_experiment/fall_2015/")

