all: libferminode.so
	#output.mp4

libferminode.so:
	python setup.py build_ext --inplace

output.mp4:
	ffmpeg -i out.mp4 -vf "setpts=20.0*PTS" output.mp4

clean:
	trash libferminode.so libferminode.c output.mp4
