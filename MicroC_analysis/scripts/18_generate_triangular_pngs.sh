size=800

for inFile in $(ls -1 *Fig5*toModify*.png 2> /dev/null);
do
    outFile=${inFile%_toModify.png}.png
    echo "Input file : ${outFile}"
    echo "Cropping..."
    convert -background white -alpha remove ${inFile} -crop ${size}x0 ${outFile}

    if [[ ! -e ${outFile%.png}-0.png ]];
    then
	for file in $(ls -1 ${outFile%.png}.png) ;
	do
	    convert +repage ${file} _tmp ; mv _tmp ${file}
	    identify ${file}
	    echo "Trimming..."
	    convert -trim ${file} ${file%.png}_trimmed.png
	    convert +repage ${file%.png}_trimmed.png _tmp ; mv _tmp ${file%.png}_trimmed.png
	    identify ${file%.png}_trimmed.png
	    echo "Rotating..."
	    convert -rotate 45 -background white -alpha remove ${file%.png}_trimmed.png _tmp_${file}
	    convert +repage _tmp_${file} _tmp; mv _tmp _tmp_${file}
	    identify _tmp_${file}
	    echo "Cropping..."
	    convert -background white -alpha remove _tmp_${file} -crop 1042x521-0-0 ${file%.png}_triang.png
	    convert +repage ${file%.png}_triang.png _tmp ; mv _tmp ${file%.png}_triang.png
	    identify ${file%.png}_triang.png
	    rm -fr ${file} _tmp*png *_trimmed.png
	done
    fi  
    
    if [[ -e ${outFile%.png}-0.png ]];
    then
	for file in $(ls -1 ${outFile%.png}-?.png) ;
	do
	    convert +repage ${file} _tmp ; mv _tmp ${file}
	    identify ${file}
	    echo "Trimming..."
	    convert -trim ${file} ${file%.png}_trimmed.png
	    convert +repage ${file%.png}_trimmed.png _tmp ; mv _tmp ${file%.png}_trimmed.png
	    identify ${file%.png}_trimmed.png
	    echo "Rotating..."
	    convert -rotate 45 -background white -alpha remove ${file%.png}_trimmed.png _tmp_${file}
	    convert +repage _tmp_${file} _tmp; mv _tmp _tmp_${file}
	    identify _tmp_${file}
	    echo "Cropping..."
	    convert -background white -alpha remove _tmp_${file} -crop 1042x521-0-0 ${file%.png}_triang.png
	    convert +repage ${file%.png}_triang.png _tmp ; mv _tmp ${file%.png}_triang.png
	    identify ${file%.png}_triang.png
	    rm -fr ${file} _tmp*png *_trimmed.png
	done
    fi

    rm -fr ${outFile} #${inFile}
done


for annotation in FALSE TRUE;
do
    n=0
    for region in chr2L_16000000bp_17000000bp chr2L_16321000bp_16512000bp chr2R_11000000bp_12000000bp chr2R_11421500bp_11612500bp chrX_14600000bp_14800000bp chr3R_16625000bp_17000000bp ;
    do
	condition=Embryo_WTGAGA
	mv -v scoreMapsk250kexp500_${region}_r1000bp_Fig2_annot_${annotation}-0_triang.png scoreMapsk250kexp500_${region}_r1000bp_Fig2_annot_${annotation}_${condition}_triang.png
	condition=Embryo_GAGA14
	mv -v scoreMapsk250kexp500_${region}_r1000bp_Fig2_annot_${annotation}-1_triang.png scoreMapsk250kexp500_${region}_r1000bp_Fig2_annot_${annotation}_${condition}_triang.png
	condition=Embryo_GAGA34
	mv -v scoreMapsk250kexp500_${region}_r1000bp_Fig2_annot_${annotation}-2_triang.png scoreMapsk250kexp500_${region}_r1000bp_Fig2_annot_${annotation}_${condition}_triang.png
	condition=Embryo_GAGAmut
	mv -v scoreMapsk250kexp500_${region}_r1000bp_Fig2_annot_${annotation}-3_triang.png scoreMapsk250kexp500_${region}_r1000bp_Fig2_annot_${annotation}_${condition}_triang.png
	
	condition=Embryo_WTGAGA
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig2_annot_${annotation}-0_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig2_annot_${annotation}_${condition}_triang.png
	condition=Embryo_GAGA14
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig2_annot_${annotation}-1_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig2_annot_${annotation}_${condition}_triang.png
	condition=Embryo_GAGA34
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig2_annot_${annotation}-2_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig2_annot_${annotation}_${condition}_triang.png
	condition=Embryo_GAGAmut
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig2_annot_${annotation}-3_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig2_annot_${annotation}_${condition}_triang.png
    done
done

for annotation in FALSE TRUE;
do
    n=0
    for region in chr2L_16000000bp_17000000bp chr2L_16321000bp_16512000bp chr2R_11000000bp_12000000bp chr2R_11421500bp_11612500bp chrX_14600000bp_14800000bp chr3R_16625000bp_17000000bp ;
    do
	condition=ED_PH18
	mv -v scoreMapsk250kexp500_${region}_r1000bp_Fig4_annot_${annotation}-0_triang.png scoreMapsk250kexp500_${region}_r1000bp_Fig4_annot_${annotation}_${condition}_triang.png
	condition=ED_PH29
	mv -v scoreMapsk250kexp500_${region}_r1000bp_Fig4_annot_${annotation}-1_triang.png scoreMapsk250kexp500_${region}_r1000bp_Fig4_annot_${annotation}_${condition}_triang.png

	condition=ED_PH18
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig4_annot_${annotation}-0_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig4_annot_${annotation}_${condition}_triang.png
	condition=ED_PH29
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig4_annot_${annotation}-1_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig4_annot_${annotation}_${condition}_triang.png	
    done
done

for annotation in FALSE TRUE;
do
    n=0
    for region in chr2L_16000000bp_17000000bp chr2L_16321000bp_16512000bp chr2R_11000000bp_12000000bp chr2R_11421500bp_11612500bp chrX_14600000bp_14800000bp chr3R_16625000bp_17000000bp ;
    do
	condition=larvae_DWT
	mv -v scoreMapsk250kexp500_${region}_r1000bp_Fig5_annot_${annotation}-0_triang.png scoreMapsk250kexp500_${region}_r1000bp_Fig5_annot_${annotation}_${condition}_triang.png
	condition=LD_DPRE2Fab7
	mv -v scoreMapsk250kexp500_${region}_r1000bp_Fig5_annot_${annotation}-1_triang.png scoreMapsk250kexp500_${region}_r1000bp_Fig5_annot_${annotation}_${condition}_triang.png
	condition=LD_DPRE2en
	mv -v scoreMapsk250kexp500_${region}_r1000bp_Fig5_annot_${annotation}-2_triang.png scoreMapsk250kexp500_${region}_r1000bp_Fig5_annot_${annotation}_${condition}_triang.png
        condition=LD_DPRE2hsp26Pho
	mv -v scoreMapsk250kexp500_${region}_r1000bp_Fig5_annot_${annotation}-3_triang.png scoreMapsk250kexp500_${region}_r1000bp_Fig5_annot_${annotation}_${condition}_triang.png
        condition=LD_DPRE2hsp26
	mv -v scoreMapsk250kexp500_${region}_r1000bp_Fig5_annot_${annotation}-4_triang.png scoreMapsk250kexp500_${region}_r1000bp_Fig5_annot_${annotation}_${condition}_triang.png
	condition=LD_DPRE2
	mv -v scoreMapsk250kexp500_${region}_r1000bp_Fig5_annot_${annotation}-5_triang.png scoreMapsk250kexp500_${region}_r1000bp_Fig5_annot_${annotation}_${condition}_triang.png
	condition=LD_DPRE2Virilis
	mv -v scoreMapsk250kexp500_${region}_r1000bp_Fig5_annot_${annotation}-6_triang.png scoreMapsk250kexp500_${region}_r1000bp_Fig5_annot_${annotation}_${condition}_triang.png
	condition=LD_DPRE2x3end
	mv -v scoreMapsk250kexp500_${region}_r1000bp_Fig5_annot_${annotation}-7_triang.png scoreMapsk250kexp500_${region}_r1000bp_Fig5_annot_${annotation}_${condition}_triang.png	
	condition=larvae_DWT
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}-0_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}_${condition}_triang.png
	condition=LD_DPRE2Fab7
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}-1_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}_${condition}_triang.png
	condition=LD_DPRE2en
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}-2_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}_${condition}_triang.png
        condition=LD_DPRE2hsp26Pho
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}-3_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}_${condition}_triang.png
        condition=LD_DPRE2hsp26
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}-4_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}_${condition}_triang.png
	condition=LD_DPRE2
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}-5_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}_${condition}_triang.png
	condition=LD_DPRE2Virilis
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}-6_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}_${condition}_triang.png
	condition=LD_DPRE2x3end
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}-7_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}_${condition}_triang.png	
    done
done

for annotation in FALSE TRUE;
do
    n=0
    for region in chr2L_16000000bp_17000000bp chr2L_16321000bp_16512000bp;
    do
	condition=larvae_DWT
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}-0_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}_${condition}_triang.png
	condition=LD_DPRE2Fab7
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}-1_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}_${condition}_triang.png
	condition=LD_DPRE2en
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}-2_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}_${condition}_triang.png
        condition=LD_DPRE2hsp26Pho
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}-3_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}_${condition}_triang.png
        condition=LD_DPRE2hsp26
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}-4_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}_${condition}_triang.png
	condition=LD_DPRE2
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}-5_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}_${condition}_triang.png
	condition=LD_DPRE2Virilis
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}-6_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}_${condition}_triang.png
	condition=LD_DPRE2x3end
	mv -v scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}-7_triang.png scoreMapsk250kexp500_${region}_r3000bp_Fig5_annot_${annotation}_${condition}_triang.png	
    done
done
