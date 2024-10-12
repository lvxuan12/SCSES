#/bin/bash
#bash ~/.bashrc
#source activate R_env
R -e ".libPaths('/home/Liulab/wenx/R/x86_64-pc-linux-gnu-library/4.2/');options(browser = 'firefox');shiny::runApp('/disk/lvxuan/Single-Splicing/src/app/',host='192.168.117.77',port=9999)"
