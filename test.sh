file=$HOME/DCP/betty_stereo.wav
src/leqm-nrt $file -chconfcal 0 0 > test
diff -u ref test && echo "PASS" || echo "FAIL"
