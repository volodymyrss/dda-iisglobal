from dataanalysis import analysisfactory

#analysisfactory.AnalysisFactory.get_by_name('ghost_bustersSpectra').cached=True

c = analysisfactory.AnalysisFactory.get_by_name('ii_spectra_extract')

c.__class__.fullbkg=True
#c.promote()

