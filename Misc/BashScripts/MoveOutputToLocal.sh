#!/bin/bash

LocalPath='../../Output/'

Host='s113245@login.gbar.dtu.dk'
RemotePath='~/BAM_artikel/Results/'

scp -r $Host:$RemotePath/* $LocalPath 
