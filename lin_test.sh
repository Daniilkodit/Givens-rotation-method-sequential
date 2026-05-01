#!/bin/bash
prog="./a.out"

if [[ -x ${prog} ]]; then
        for ((m = 3; m <= 200; m += 1))
                do echo "=========================== n = 2000, m = $m, s = 1  ===========================" 
                   ${prog} 2000 $m 1 1    # ← ИСПРАВЛЕНО: r = n
        done
               
fi
