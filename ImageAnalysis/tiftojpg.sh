for img in `ls *.tif`
do
  convert $img $img.jpg
done
