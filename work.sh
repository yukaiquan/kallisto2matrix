# copy from project
cp ../target/release/kallisto2matrix.exe ./
cp ../target/x86_64-unknown-linux-musl/release/kallisto2matrix ./

# run example
cd example
../kallisto2matrix -i samples.txt -o test

# run example of salmon
../kallisto2matrix.exe -i ./samples_salmon.txt -o salmon -t salmon
