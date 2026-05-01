#!/bin/bash
echo "=== Testing matrices with block size 2 ==="
echo "a.txt"; ./a.out 4 2 4 0 a.txt; echo
echo "a20.txt"; ./a.out 4 2 4 0 a20.txt; echo  
echo "b.txt"; ./a.out 4 2 4 0 b.txt; echo
echo "c.txt"; ./a.out 6 2 6 0 c.txt; echo
echo "d.txt"; ./a.out 6 2 6 0 d.txt; echo
echo "e.txt"; ./a.out 6 2 6 0 e.txt; echo
echo "f.txt"; ./a.out 4 2 4 0 f.txt; echo