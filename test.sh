file=$HOME/DCP/betty_stereo.wav
LD_LIBRARY_PATH=build/src build/src/leqm-nrt $file -chconfcal 0 0 > test
diff -u ref test && echo "PASS" || echo "FAIL"
