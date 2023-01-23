#!/bin/bash

nomigeneral=$1

mencoder "mf://$nomigeneral*png" -mf fps=6 -o prueba.avi -ovc lavc -lavcopts vcodec=mpeg4:vbitrate=800
