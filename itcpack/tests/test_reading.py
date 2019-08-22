from itcpack import readitc

def test_readitc_onvpitc_headerfields():
    '''header fields in .itc file contain cell volume, starting cell concentration
    and syringe concentration
    '''
    itcdata=readitc.readitc("datafiles/040419b.itc")

    #syringe conc in M and vo (cell volume) in L pulled from data file
    assert itcdata.syrconc==.00202
    assert itcdata.vo==0.0014207

    #additionally, starting mconc in cell should be stored in dataclass field expdetails
    known_expdetails=readitc.expdetails(.00202,.00004,.0014207)
    assert itcdata.expdetails==known_expdetails

def test_readitc_onvpitc_injections_extended_data():
    '''test that list of injections is read in properly'''
    itcdata=readitc.readitc("datafiles/040419b.itc")

    #there are 21 injections in this file + pre-injection regime. Test size and injnum labels
    assert len(itcdata.injections)==22
    assert [x.injnum for x in itcdata.injections]==list(range(22))

    #test that lengths are correct for injection at start, middle,end
    #start with pre-injection baseline ('faux injection')
    fauxinj=itcdata.injections[0]
    assert len(fauxinj.seconds_total)==150 and len(fauxinj.power)==150 and len(fauxinj.seconds)==150
    
    firstinj=itcdata.injections[1]
    assert len(firstinj.seconds_total)==120 and len(firstinj.power)==120 and len(firstinj.seconds)==120
    
    midinj=itcdata.injections[15]
    assert len(midinj.seconds_total)==105 and len(midinj.power)==105 and len(midinj.seconds)==105
    
    lastinj=itcdata.injections[21]
    assert len(lastinj.seconds_total)==79 and len(lastinj.power)==79 and len(lastinj.seconds)==79

    #make sure seconds tally and seconds_total are correct
    cumulative_seconds=0
    for injidx in range(len(itcdata.injections)):
        cumulative_seconds+= itcdata.injections[injidx].seconds_total[0]-cumulative_seconds
        assert cumulative_seconds==itcdata.injections[injidx].seconds_total[0]
        assert itcdata.injections[injidx].seconds_total[-1]-itcdata.injections[injidx].seconds_total[0]\
                ==itcdata.injections[injidx].seconds[-1]
        cumulative_seconds+=itcdata.injections[injidx].seconds_total[-1]

def test_readitc_onvpitc_injections_minimal_data():
    '''test that list of injections is read in properly. This file only collected 3 values throughout titration'''
    itcdata=readitc.readitc("datafiles/061810b.itc")

    #there are 24 injections in this file + pre-injection range Test size and injnum labels
    assert len(itcdata.injections)==25
    assert [x.injnum for x in itcdata.injections]==list(range(25))

    #test that lengths are correct for injection at start, middle,end
    #start with pre-injection baseline ('faux injection')
#    fauxinj=itcdata.injections[0]
#    assert len(fauxinj.seconds_total)==150 and len(fauxinj.power)==150 and len(fauxinj.seconds)==150
#    
#    firstinj=itcdata.injections[1]
#    assert len(firstinj.seconds_total)==120 and len(firstinj.power)==120 and len(firstinj.seconds)==120
#    
#    midinj=itcdata.injections[15]
#    assert len(midinj.seconds_total)==105 and len(midinj.power)==105 and len(midinj.seconds)==105
#    
#    lastinj=itcdata.injections[21]
#    assert len(lastinj.seconds_total)==79 and len(lastinj.power)==79 and len(lastinj.seconds)==79

    #make sure seconds tally and seconds_total are correct
    cumulative_seconds=0
    for injidx in range(len(itcdata.injections)):
        cumulative_seconds+= itcdata.injections[injidx].seconds_total[0]-cumulative_seconds
        assert cumulative_seconds==itcdata.injections[injidx].seconds_total[0]
        assert itcdata.injections[injidx].seconds_total[-1]-itcdata.injections[injidx].seconds_total[0]\
                ==itcdata.injections[injidx].seconds[-1]
        cumulative_seconds+=itcdata.injections[injidx].seconds_total[-1]


    #add test for mtoti,mtotf,ltoti,ltotf ?

#def test_create_titration_dataset():
#    itcdata=readitc.readitc("datafiles/040419b.itc")
#    itcdata.create_titration_dataset()
#
#
#    #syringe conc in M and vo (cell volume) in L pulled from data file
#    assert itcdata.syrconc==.00202
#    assert itcdata.vo==0.0014207

