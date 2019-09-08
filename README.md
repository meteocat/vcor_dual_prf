{\rtf1\ansi\ansicpg1252\cocoartf1504\cocoasubrtf830
{\fonttbl\f0\fmodern\fcharset0 Courier;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
\paperw11900\paperh16840\margl1440\margr1440\vieww10700\viewh8400\viewkind0
\deftab709
\pard\pardeftab709\sl302\pardirnatural\qc\partightenfactor0

\f0\fs24 \cf0 Correction of dual-PRF velocity dealiasing errors\cf2 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 \
\pard\pardeftab709\sl282\sa170\pardirnatural\qc\partightenfactor0
\cf2 =====================================================================\
\pard\pardeftab709\sl311\sa170\pardirnatural\qj\partightenfactor0

\fs22 \cf0 \kerning1\expnd0\expndtw0 \outl0\strokewidth0 The dual-PRF dealiasing errors constitute a serious handicap for the global dealiasing of the radial velocity field and for derived applications such as wind field estimation and meteorological signature detection. \
This repository includes a correction function openly available for the weather radar community. The function allows the user to tailor the dual-PRF error correction by specifying the neighbour kernel and selecting the statistic used for the estimation of the reference velocity. The correction procedures proposed in the literature may be reproduced through the particular selection of these statistics: the mean velocity (Joe and May, 2003), the median velocity (Holleman and Beekhuis, 2003), the circular mean of the PRF-scaled phase (Altube et al., 2017) and the circular mean of the phase (Hengstebeck et al., 2018).}