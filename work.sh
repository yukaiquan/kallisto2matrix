# copy from project
cp ../target/release/kallisto2matrix.exe ./
cp ../target/x86_64-unknown-linux-musl/release/kallisto2matrix ./

# run example
cd example
../kallisto2matrix -i samples.txt -o test

