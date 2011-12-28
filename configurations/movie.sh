for f in *.gray; do convert -size 512x512 -depth 16 -define quantum:format=unsigned-integer -endian LSB -quality 100 -level 3800,8000 $f ${f/gray/tga}; done

convert -size 512x512 -depth 16 -define quantum:format=unsigned-integer -endian LSB -quality 100 -level 3800,8000 image.a.?.gray image.a.??.gray image.a.???.gray -depth 8 -type scaled img%03d.sgi

animate -delay 10 img???.png

mencoder "mf://*.jpg" -mf fps=10 -o test.avi -ovc lavc -lavcopts vcodec=msmpeg4v2:vbitrate=800

mencoder "mf://*.?.jpg" "mf://*.??.jpg" "mf://*.???.jpg" -quiet -mf fps=30 -o chem.avi -ovc lavc -lavcopts vcodec=msmpeg4v2

mencoder "mf://image.b.?.gray.jpg" "mf://image.b.??.gray.jpg" -quiet -o testmp4.avi -ovc x264 -x264encopts threads=auto:bitrate=10000

#Concatenate several movie files into one!!!
mencoder *.avi -o final.avi -ovc copy
 
mplayer -loop 0 -vf dsize=3100:-2 -geometry 0:0 test.avi

FROM: http://www.tevs.eu/blog_8_comments.html
Encodes png sequence into quicktime compatible mpeg4
mencoder mf://PATH/*.png -mf fps=FPS -vf scale=720x576,harddup -ovc x264 -x264encopts pass=1:bitrate=BITRATE:threads=auto:frameref=5:bframes=0:nob_pyramid:nob_adapt:direct_pred=auto:subq=6:mixed_refs:nodct_decimate:no_psnr -nosound -ofps FPS -noskip -of rawvideo -o /dev/null
mencoder mf://PATH/*.png -mf fps=FPS -vf scale=720x576,harddup -ovc x264 -x264encopts pass=2:bitrate=BITRATE:threads=auto:frameref=5:bframes=0:nob_pyramid:nob_adapt:direct_pred=auto:subq=6:mixed_refs:nodct_decimate:no_psnr -nosound -ofps FPS -noskip -of rawvideo -o output.264

FROM: http://wiki.qwdrama.com/Mencoder_howto
xVid
mencoder "mf://*.tga" -mf fps=25 -o /dev/null -ovc xvid -xvidencopts pass=1:bitrate=3000 -vf scale=320:240,eq2=1.5
mencoder "mf://*.tga" -mf fps=25 -o output.avi -ovc xvid -xvidencopts pass=2:bitrate=3000 -vf scale=320:240,eq2=1.5

FROM: http://www.entropyfarm.org/software/tricks/
converts a movie into individual frames, 'n' is from 1-9 for compression level
mplayer -nosound -vo png:z=n inputfile

FROM: http://personal.cscs.ch/~mvalle/mencoder/mencoder.html
The 50 factor can vary between 40 and 60 to trade quality for size.
optimal_bitrate = 50 * 25 * width * height / 256
For the mpeg4 codec:
opt="vbitrate=2160000:mbd=2:keyint=132:v4mv:vqmin=3:vlelim=-4:vcelim=7:lumi_mask=0.07:dark_mask=0.10:naq:vqcomp=0.7:vqblur=0.2:mpeg_quant"
For scientific visualization movies (lot of hard edges) I have found interesting also the following set of parameters:
opt="vbitrate=2160000:mbd=2:keyint=132:v4mv:vqmin=3:lumi_mask=0.07:dark_mask=0.2:mpeg_quant:scplx_mask=0.1:tcplx_mask=0.1:naq"
Anyway, I plan to test various s of parameters to find the best one at least on common (for me) cases .
For the msmpeg4v2 codec:
opt="vbitrate=2160000:mbd=2:keyint=132:vqblur=1.0:cmp=2:subcmp=2:dia=2:mv0:last_pred=3"


FROM: http://www4.mplayerhq.hu/DOCS/HTML/en/menc-feat-enc-libavcodec.html
Very high quality	vcodec=mpeg4:mbd=2:mv0:trell:v4mv:cbp:last_pred=3:predia=2:dia=2:vmax_b_frames=2:vb_strategy=1:precmp=2:cmp=2:subcmp=2:preme=2:qns=2	6fps	0dB
High quality	vcodec=mpeg4:mbd=2:trell:v4mv:last_pred=2:dia=-1:vmax_b_frames=2:vb_strategy=1:cmp=3:subcmp=3:precmp=0:vqcomp=0.6:turbo	15fps	-0.5dB
Fast	vcodec=mpeg4:mbd=2:trell:v4mv:turbo	42fps	-0.74dB
Realtime	vcodec=mpeg4:mbd=2:turbo

trell - Trellis searched quantization.  This will find the optimal encoding for each 8x8 block. 
mbd=2 - Macroblock  decision  algorithm  (high quality mode), encode each macro block in all modes and choose the best. ('2' - Select the MB mode which has the best rate distortion.)
v4mv - Allow 4 motion vectors per macroblock (slightly better quality).  Works better if used with mbd>0.
dia=2 - Diamond type & size for motion estimation (-99-6). Bigger diamonds allow a wider search for the best motion vector
keyint=30 - maximum interval between keyframes in frames (0-300)


v4mv:mbd=2:trell:dia=2


http://www.ffmpegx.com/options.html

FROM: http://lists.mplayerhq.hu/pipermail/mencoder-users/2005-April/000783.html
qns=3:cmp=2:subcmp=6:preme=2:predia=3:dia=3:last_pred=3:vmax_b_frames=1:mv0:cbp:qprd:trell:preme=2:keyint=300:v4mv:mbd=2:vqcomp=0.7
FROM: http://lists.mplayerhq.hu/pipermail/mencoder-users/2005-April/000791.html
-lavcopts vcodec=mpeg4:mbd=2:v4mv:trell:vpass=$N:vbitrate=1800 -ofps 24000/1001
FROM: http://lists.mplayerhq.hu/pipermail/mencoder-users/2005-April/000794.html
vcodec=mpeg4:vbitrate=XXX:vqmin=1:lmin=1:v4mv:mbd=2:trell:subcmp=2:cmp=2:precmp=2:mv0:autoaspect:vpass=1:mpeg_quant:psnr:vqcomp=0.7:vqblur=0.3

convert -size 256x256 -depth 16 -define quantum:format=unsigned-integer -endian LSB -quality 100 -level 1442,12815 image.a.?????.gray -resize 300% img%03d.jpg

FPS=10
OUTPUT=hello

mencoder "mf://*.jpg" -mf fps=$FPS -ovc lavc -lavcopts vcodec=msmpeg4v2:keyint=1 -o $OUTPUT.avi

