xhost + 192.168.1.15
docker run -it -e DISPLAY=192.168.1.15:0 -v /Users/ach3377:/home/neuro achennings/afni_fsl_py 