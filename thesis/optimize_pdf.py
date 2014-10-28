#!/usr/bin/python

import os

for fich in os.listdir('pdf/'):
    if fich[-3:]=="pdf":
        print(fich)
        os.system("gs -sDEVICE=pdfwrite -dCompatibilityLevel=1.4 -dPDFSETTINGS=/printer -dUseCIEColor -dNOPAUSE -dQUIET -dBATCH -sOutputFile=pix/%s pdf/%s" % (fich,fich))
