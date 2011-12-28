#!/bin/bash

grep --color=always -H -n -C 5 ERROR *|less -R
