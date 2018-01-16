./clean.sh
./compile.sh
./test.sh | grep -v VERBOSE > toto
vimdiff toto expected_output
rm -f toto
