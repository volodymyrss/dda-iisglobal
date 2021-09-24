import os
import glob
import pilton

from ddosa import (DataAnalysis, 
                  ghost_bustersSpectra, 
                  CatForSpectraFromImaging, 
                  IBIS_ICRoot,
                  ISGRIResponse,
                  ScWData,
                  BinMapsSpectra,
                  ibis_gti,
                  construct_gnrl_scwg_grp,
                  import_attr,
                  set_attr,
                  DataFile,
                  remove_withtemplate,
                  heatool
                  )


class ii_spectra_extract(DataAnalysis):
    input_gb=ghost_bustersSpectra
    input_cat=CatForSpectraFromImaging # or virtual
    input_ic=IBIS_ICRoot
    input_response=ISGRIResponse
    input_scw=ScWData()
    input_maps=BinMapsSpectra

    input_gti=ibis_gti

    cached=True

    version="v3"

    #input_bins=SpectraBins
    #input_cat=CatExtract
    #input_imgconfig=ImagingConfig

    shdtype="BIN_S"
    binary="ii_spectra_extract"

    fullbkg=False

    usebkg=True

    copy_cached_input=False

    def get_version(self):
        v =  super().get_version()
        if self.fullbkg:
            v += ".fullbkg"
        return v

    def morphed_gb(self):
        return self.input_gb.corshad.path

        f = fits.open(self.input_gb.corshad.path)
        fn = "morphed_gb.fits"

        crab_pif = fits.open('/home/savchenk/work/integral/spectral_background/isgri_pif.fits')[2].data
        bkg = fits.open(self.morph_bkg())

        base = None
        for e in f[2:]:
            if e.header['E_MIN'] == 32.0:
                print('\033[31mBASE: \033[0m', e.header['E_MIN'], e.header['E_MAX'], np.quantile(e.data.flatten(), 0.99))
                if e.header['SHD_TYPE'] == 'EFFICIENCY':
                    # e.data[e.data > 0] = np.ones_like(e.data)[e.data > 0]
                    e.data = 1 * np.ones_like(e.data)
                elif e.header['SHD_TYPE'] == 'DETECTOR':
                    # base = e.data.copy()
                    # e.data = 100 * np.ones_like(e.data)

                    e.data = 50 * crab_pif.copy() + 1000 * bkg[2].data # all the same for now
                    # e.data = 7 * bkg[30].data # all the same for now
                elif e.header['SHD_TYPE'] == 'VARIANCE':
                    # base = e.data.copy()
                    # e.data = 100 * np.ones_like(e.data)

                    e.data = 0 * crab_pif.copy() + 7 * bkg[2].data # all the same for now
                    # e.data = 7 * bkg[30].data # all the same for now

            # elif e.header['E_MIN'] == 32.0:
            #     print('\033[31mTEST: \033[0m', e.header['E_MIN'], e.header['E_MAX'], np.quantile(e.data.flatten(), 0.99))
            #     if e.header['SHD_TYPE'] == 'EFFICIENCY':
            #         #e.data = 0.5 * np.ones_like(e.data)
            #         e.data = 1 * np.ones_like(e.data)
            #     elif e.header['SHD_TYPE'] == 'DETECTOR':
            #         # e.data = base.copy()
            #         # e.data = 50 * np.ones_like(e.data)
            #         e.data = 20 * crab_pif.copy() + 7 * bkg[30].data
            #         # e.data = 5 * bkg[30].data
            #     elif e.header['SHD_TYPE'] == 'VARIANCE':
            #         # e.data = base.copy()
            #         # e.data = 50 * np.ones_like(e.data)
            #         e.data = 20 * crab_pif.copy() + 7 * bkg[30].data
            #         # e.data = 5 * bkg[30].data
            else:
                e.data = np.zeros_like(e.data)

        f.writeto(fn, overwrite=True)
        return fn

    def morph_bkg(self):
        p = self.input_maps.back.get_path()
        return p

        f = fits.open(p)

        print("morph_bkg:", p)

        for e in f[2:]:
            s1, s2 = e.data.shape 
            g1, g2 = np.meshgrid(np.arange(s2), np.arange(s1))
            #e.data[:, :] = 10 * (1 - ((g1[:, :] - 134/2)/(134/2))**4 * ((g2[:, :] - 130/2)/(130/2))**4)
            e.data[:, :] = g1[:, :]
            e.data[:, :] /= e.data[:, :].mean() * 3
            
            # e.data[:, :] = 55* np.ones_like(g1)
            print("G1:", g1)

        fn = "morphed_bkg.fits"
        f.writeto(fn, overwrite=True)

        return fn

    def morph_corr(self):
        p = self.input_maps.corr.get_path()
        return p

        f = fits.open(p)

        print("morph_corr:", p)

        for e in f[2:]:
            e.data[:, :] = np.ones_like(e.data)
            
        fn = "morphed_corr.fits"
        f.writeto(fn, overwrite=True)

        return fn

    def main(self):
        # if self.fullbkg:
        #     pass
        #     #self.binary = "/home/savchenk/work/integral/spectral_background/ii_spectra_extract/ii_spectra_extract"

        if hasattr(self.input_cat,'empty_results') and self.input_cat.empty_results:
            print("empty here")
            self.empty_results=True
            return

        att=self.input_scw.auxadppath+"/attitude_historic.fits"
        if os.path.exists(att):
            att=self.input_scw.auxadppath+"/attitude_historic.fits[AUXL-ATTI-HIS,1,BINTABLE]"
            attp=att
        else:
            att=self.input_scw.auxadppath+"/attitude_snapshot.fits[AUXL-ATTI-SNA,1,BINTABLE]"
            attp_fn=glob.glob(self.input_scw.auxadppath+"/attitude_predicted_*.fits*")[0]
            attp=attp_fn+"[AUXL-ATTI-PRE,1,BINTABLE]"
        
        construct_gnrl_scwg_grp(self.input_scw,[\
                    #self.input_gb.corshad.path,
                    self.morphed_gb(),
                    self.input_scw.auxadppath+"/time_correlation.fits[AUXL-TCOR-HIS]",
                    att,
                    attp,
                    self.input_gti.output_gti.path
                ])
                    #self.input_cat.cat.path,

        import_attr(self.input_scw.scwpath+"/swg.fits",["OBTSTART","OBTEND","TSTART","TSTOP","SW_TYPE","TELAPSE","SWID"])
        set_attr({'ISDCLEVL':"BIN_S"})
        #set_attr({'INSTRUME':"IBIS"},"og.fits")

        #construct_gnrl_scwg_grp_idx(self.input_scw,[\
        #            "og.fits",
        #        ])
        #set_attr({'ISDCLEVL':"BIN_I"},"og_idx.fits")

      #  construct_og(self.input_scw,[\
      #              "og_idx.fits",
      #          ])
      #  set_attr({'ISDCLEVL':"BIN_I"},"ogg.fits")

        #remove_withtemplate("isgri_srcl_res.fits(ISGR-SRCL-RES.tpl)")

        pif_fn,pif_tpl="isgri_pif.fits","(ISGR-PIF.-SHD-IDX.tpl)"
        spec_fn,spec_tpl="isgri_spectrum.fits","(ISGR-EVTS-SPE-IDX.tpl)"

        remove_withtemplate(pif_fn+pif_tpl)
        remove_withtemplate(spec_fn+spec_tpl)

        ht=heatool(self.binary)
        ht['outSwg']="og.fits"
        ht['inCat']=self.input_cat.cat.get_path()
        ht['outPif']=pif_fn+pif_tpl
        ht['outSpec']=spec_fn+spec_tpl
        ht['mask']=self.input_ic.ibisicroot+"/mod/isgr_mask_mod_0003.fits[ISGR-MASK-MOD,1,IMAGE]"
        ht['tungAtt']=self.input_ic.ibisicroot+"/mod/isgr_attn_mod_0010.fits[ISGR-ATTN-MOD,1,BINTABLE]"
        ht['aluAtt']=self.input_ic.ibisicroot+"/mod/isgr_attn_mod_0011.fits[ISGR-ATTN-MOD,1,BINTABLE]"
        ht['leadAtt']=self.input_ic.ibisicroot+"/mod/isgr_attn_mod_0012.fits[ISGR-ATTN-MOD,1,BINTABLE]"
        ht['idx_isgrResp']=self.input_response.rmf_path
        ht['isgrUnifDol']=self.input_maps.unif.path
        if self.usebkg:
            # ht['isgrBkgDol']=self.input_maps.back.path
            ht['isgrBkgDol']=self.morph_bkg()
        else:
            ht['isgrBkgDol']="-"
        ht['corrDol']=self.morph_corr()
        ht['OutType']=self.shdtype
        if self.fullbkg:
            ht['method_cor'] = 4
        else:
            ht['method_cor'] = 1 # global correction

        if hasattr(self,'input_bins') and not self.input_bins.rmfbins:
            ebins=self.input_bins.bins
            ht['num_band']=len(ebins)
            ht['E_band_min'],ht['E_band_max']=[" ".join(["%.5lg"%b for b in a]) for a in zip(*ebins)]
        else:
            rmf=self.input_bins.get_binrmfext() if hasattr(self,'input_bins') else self.input_response.rmf_path
            ht['num_band'] = -1
            ht['idx_isgrResp'] = rmf #  +"[1]" will this work?


        #for k in ['SearchMode','ToSearch','CleanMode','MinCatSouSnr','MinNewSouSnr','NegModels','DoPart2']: # dopart2 is flow control, separately
        #    ht[k]=getattr(self.input_imgconfig,k)

        try:
            ht.run()
            self.spectrum=DataFile(spec_fn)
            self.pifs=DataFile(pif_fn)
        except pilton.HEAToolException as e:
            self.empty_results=True


