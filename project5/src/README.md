To compile and link the files you use the `makefile` by typing e.g.

```
make -f makefile SIRS_test
```

The program writes to standard output the temporal estimates of S, I and R. To run the file for b=1, and re-direct the results to a file called `sir_test_b=1.out`, type

```
./SIRS_test 1.0  > sirs_b=1.out

```

