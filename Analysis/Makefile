CC = clang
CFLAGS = -g -DPPOS
LDFLAGS = -g

PYTHON = python

#Make a movie out of the snapshots using ffmpeg
MOVIENAME = movie.mp4
FRAMEDIR = frames
FRAMEBASE = frame
FRAMESPERSECOND = 6

read_snapshot: read_snapshot.o

read_snapshot.o: read_snapshot.c

.PHONY: clean
clean:
	rm -rf *.o *.dSYM *.mp4 read_snapshot

.PHONY: movie
movie: pngframes 
	ffmpeg -f image2 -r $(FRAMESPERSECOND) -i $(FRAMEDIR)/$(FRAMEBASE)%03d.png -vcodec mpeg4 -y $(MOVIENAME)

pngframes: read_snapshot read_snapshot.sh 3dsnapshot.py 
	rm -rf pngframes
	$(PYTHON) 3dsnapshot.py $(FRAMEDIR) $(FRAMEBASE)
	touch pngframes