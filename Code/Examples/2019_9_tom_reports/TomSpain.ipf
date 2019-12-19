#pragma rtGlobals=3		// Use modern global access method and strict wave access.
#pragma ModuleName = ModTomSpain

#include  ":::Lib:UtilIgorPro:Util:PlotUtil"
Static function plot_single(wave_x,wave_y,color)
	Wave wave_x, wave_y
	String color
	ModPlotUtil#PlotWithFiltered(wave_y,X=wave_x,nFilterPoints=50,color=color,alpha=0.9,plot_raw=1)
	DoUpdate
End Function

Static Function /Wave duplicated(wave_tmp)
	Wave wave_tmp
	String NameProc = (NameOfWave(wave_tmp) + "_Pro")
	KillWaves /Z $NameProc
	Duplicate /O wave_tmp $NameProc
	Wave processed = $NameProc
	return processed
End Function

Static function /Wave process_single_y(wave_y)
	Wave wave_y
	Wave processed =  duplicated(wave_y)
	processed[] = processed[p] - 37e-12
	processed[] = -1e12 * processed[p] 
	return processed
End Function

Static function /Wave process_single_t(wave_y)
	Wave wave_y
	Wave processed =  duplicated(wave_y)
	Variable min_v = processed[0]
	processed[] = processed[p] - min_v
	processed[] = processed[p] - 0.15
	return processed
End Function


Static function /Wave process_single_x(wave_y)
	Wave wave_y
	Wave processed =  duplicated(wave_y)
	Variable min_v = processed[0]
	processed[] = processed[p] - min_v
	processed[] = processed[p] - 5e-9
	processed[] = 1e9 * processed[p]
	return processed
End Function

Static Function Main()
	ModPlotUtil#clf()
	String fig = ModPlotUtil#Figure(hide=0,heightIn=6,widthIn=14)
	Make /O/T text_wave= {"Image0863","Image0864","Image0869"}
	Make /O/T colors= {"r","b","g"}	
	Make /O offsets = {-20,5,25}
	Variable i
	String color
	String ax1 = ModPlotUtil#subplot(1,2,1)
	for (i=0; i< 3; i+= 1)
		wave t_raw = $((text_wave[i] +"Time_Ret"))
		wave F_raw = $((text_wave[i] +"Force_Ret"))
		color = colors[i]
		Wave t = process_single_t(t_raw)
		Wave F =process_single_y(F_raw)
		F[] = F[p] + offsets[i]
		plot_single(t,F,color)
	endfor
	Variable ymin = -40
	Variable ymax = 175
	ModPlotUtil#xlabel("Time (s)")
	ModPlotUtil#ylabel("Force (pN)")
	ModPlotUtil#Xlim(-0.1,0.75)
	ModPlotUtil#Ylim(ymin,ymax)
	String ax2 = ModPlotUtil#subplot(1,2,2)
	for (i=0; i< 3; i+= 1)
		wave F_raw = $((text_wave[i] +"Force_Ret"))
		wave x_raw = $((text_wave[i] +"Sep_Ret"))
		color = colors[i]
		Wave x = process_single_x(x_raw)
		Wave F =process_single_y(F_raw)
		F[] = F[p] + offsets[i]
		plot_single(x,F,color)
	endfor
	ModPlotUtil#xlabel("Extension (nm)")
	ModPlotUtil#ylabel(" ")
	ModPlotUtil#Xlim(-10,150)
	ModPlotUtil#Ylim(ymin,ymax)
	ModifyGraph noLabel(left)=1
	ModifyGraph margin(left)=-1
	DoUpdate
End Function  