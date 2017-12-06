./clean.sh
./compile.sh
./test.sh > toto
vimdiff toto expected_output
rm -f toto