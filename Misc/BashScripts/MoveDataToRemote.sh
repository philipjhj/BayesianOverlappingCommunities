#!/bin/bash

LocalPath='../../Data/Mouse*'

Host='s113245@login.gbar.dtu.dk'
RemotePath='~/BAM_artikel/Data/'

scp -r $LocalPath* $Host:$RemotePath
