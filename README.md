# ROCA
Composite Generation for GMP

## File Contents
1. ROCA.c -
   This is the main body of the code itself.
2. ROCA -
   This is the complied version of ROCA.c, feel free to compile your own after audit!
3. ROCA.log -
   The resutls of running ./ROCA will be stored here.
4. ROCAtest.log -
   The resutls of running `./ROCA -c 20` are stored here for testing.

## Compiling
ROCA.c is complied with ```gcc -I/usr/local/include ROCA.c -o ROCA -lgmp -L/usr/local/lib -std=c99```. <br>
It's dependencies are GMP 6.1.2

## Flags/Arguments
There are 7 flags for use in ROCA:
*  -b :    0 < target bit-size <= 8192  default: 512
*  -f :   fudge factor  default: 1
*  -p :   0 < prime count <= 200 default: 70
*  -c :   log2 of number of trials default: 33
*  -s :   random seed in hex default: 0x1337
*  -v :   record this number of passes default: 1
*  -o :   offset on starting k value, where offset = o*c default: 0

ROCA has some input validation on flags, and many are left as default.

## Parallelisation
There are multiple ways to parallelise this code, the offset parameter `o` can be made use of with GNU parallel
(this is what we have been doing in stress tests) as follows:

The code takes some starting value k and generates some x = kM (where M is the product of the first p primes from primecount),
such that (2x+1) and (4x+1) are of bit size `b`.

The offset paramater `o` simply shifts the starting value of k to k+offset&ast;2^{c}.

#### Example:
Say c = 20, then:<br>
`./ROCA -c 20 -o 0 ` will test the first 2^{20} k values<br>
`./ROCA -c 20 -o 1 ` will test the next 2^{20} k values<br>

#### GNU Parallel Example:
By creating a file `args.txt`:

```
-o 0
-o 1
-o 2
-o 3
-o 4
-o 5
-o 6
-o 7
-o 8
-o 9
-o 10
-o 11
```
We run the command `cat args.txt | parallel -j12 ./ROCA` to have 12 cores running 2^20 candidates each
i.e 12 &ast; 2^{20} total tests.

## Prime and Prejudice Computations
Would like to run the computation:
1. `./ROCA -b 512 -c 43 ` - This aims to find the 1024 bit example of a composite that passes 15 rounds of GMP's primality test.

2. `./ROCA -b 1024 -c 45 `- This aims to find the 2048 bit example of a composite that passes 15 rounds of GMP's primality test.

### Output

We are interested in the resulting ROCA.log file (we are looking for an entry that reads 15: x : n for some x, n).

## Test files
To check the code runs on your machine correctly, we have included the log file of running `./ROCA -c 20` for comparison. Remember to save and then remove the ROCA.log file after each computation, as the log is continually added to.
