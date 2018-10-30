### Png to mp4 ###

ffmpeg -framerate 5 -i seemsWorking%03d.png -c:v libx264 -r 30 -pix_fmt yuv420p out.mp4