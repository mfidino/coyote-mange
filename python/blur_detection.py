# ----------------------------------------------- #
# Notes:
# You'll need pip to download the necessary modules
# pip install imutils
# pip install numpy 
# pip install opencv-python

# Example of a console call 
# python blur_detection.py --images {IMAGE_PATH} --output {OUTPUTPATH.csv}
# ----------------------------------------------- #

from imutils import paths
import argparse
import cv2
import csv
 
def variance_of_laplacian(image):
  # compute the Laplacian of the image and then return the focus
  # measure, which is simply the variance of the Laplacian
  return cv2.Laplacian(image, cv2.CV_64F).var()
 
# construct the argument parse and parse the arguments
ap = argparse.ArgumentParser()
ap.add_argument("-i", "--images", required=True,
  help="path to input directory of images")
ap.add_argument("-o", "--output", required=True,
  help="path to output directory of csv file")
args = vars(ap.parse_args())


# loop over the input images
for imagePath in paths.list_images(args["images"]):
  # load the image, convert it to grayscale, and compute the
  # focus measure of the image using the Variance of Laplacian
  # method
  image = cv2.imread(imagePath)
  gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
  fm = variance_of_laplacian(gray)

  results = {'path':imagePath,'blur':fm}
  with open(args["output"],'a') as fp:
    fieldnames = ['path', 'blur']
    wr = csv.DictWriter(fp, fieldnames=fieldnames)
    fp.seek(0, 2)

    if fp.tell() == 0:
      wr.writeheader()

    wr.writerow(results)
