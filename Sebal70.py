# -*- coding: utf-8

import numpy  
import math  
import os  
import grass.script as grass
from os import system
from grass.script import core as g
  
system('clear')

MTLfile = [i for i in os.listdir('.') if i.endswith('MTL.txt')]

files = [i for i in os.listdir('.') if i.endswith('.TIF')]
basename = files[4]

u_2m = float(input('Place the wind speed value for height of the 2 m measured in the weather station - (m/s): '))

EToi = float(input('Place the instantaneous value of reference evapotranspiration (EToi) from the weather station, for the time of the satellite overpass - (mm): '))

ETo = float(input('Place the daily value of reference evapotranspiration (ETo) from the weather station - (mm): '))

g.parse_command('g.region',flags='p',rast='MDT_Sebal@PERMANENT',quiet=True)

runCC = g.parse_command('g.list',type='raster', pattern='CC_652')
runRLo = g.parse_command('g.list',type='raster', pattern='RLo')

if runCC == {}:
        print 'Importing Landsat 8 images, be patient...'
        for i in range(len(files)):
                g.parse_command('r.in.gdal', input=files[i], output=os.path.splitext(files[i])[0], overwrite=True)
        print 'Done!'
 
        print 'Calculates top-of-atmosphere reflectance and temperature for Landsat 8, be patient...'
        g.parse_command('i.landsat.toar', input=basename[0:(len(basename)-5)], output= 'LS8_corre', metfile=MTLfile, sensor='oli8',overwrite=True)
        print 'Done!'

        print 'Composite R=B6 G=B5 B=B2 Landsat 8, be patient...'
        grass.run_command('i.colors.enhance', red='LS8_corre6',green='LS8_corre5',blue='LS8_corre2',quiet=True)
        grass.run_command('r.composite',red='LS8_corre6', green='LS8_corre5', blue='LS8_corre2', output='CC_652',quiet=True,overwrite=True)
        print 'Done!'
if runRLo == {}:
        print 'Calculating NDVI...'
        grass.mapcalc('NDVI=($LS8_corre5-$LS8_corre4)/($LS8_corre5+$LS8_corre4)',
                LS8_corre5='LS8_corre5',
                LS8_corre4='LS8_corre4',
                overwrite='true',
                quiet='true')
        print 'Done!'
        
        print 'Calculating SAVI with L equal to 0.5...'
        Lsavi = 0.5 #float(input('Entre com o valor de L: '))
        grass.mapcalc('SAVI= (($LS8_corre5-$LS8_corre4)/($LS8_corre5+$LS8_corre4+$Lsavi))*(1+$Lsavi)',
                LS8_corre4='LS8_corre4',
                LS8_corre5='LS8_corre5',
                Lsavi=Lsavi,
                overwrite='true',
                quiet='true')
        print 'Done!'
        
        print 'Calculating leaf area index (LAI)...'
        grass.mapcalc('LAI=if(SAVI < 0.1, 0.00001,(if(0.1 < SAVI && SAVI < 0.687,-log((0.69-SAVI)/0.59)/0.91,if(SAVI > 0.687,6,0))))',
                overwrite='true',
                quiet='true')
        print 'Done!'
        
        print 'Calculating narrow band surface emissivity (eNBf)...'
        grass.mapcalc('eNB=if(LAI<3 && NDVI>0.,0.97+0.0033*LAI,if(LAI>=3 && NDVI>0.,0.98))',
                overwrite='true',
                quiet='true')
        grass.mapcalc('eNBf=if(eNB == 0, 0.99, eNB)',
                overwrite='true',
                quiet='true')
        print 'Done!'
        
        print 'Calculating broad band surface emissivity (e0f)...'
        grass.mapcalc('e0=if(LAI<3 && NDVI>0.,0.95+0.01*LAI,if(LAI>=3 && NDVI>0.,0.98))',
                overwrite='true',
                quiet='true') 
        grass.mapcalc('e0f=if(e0 == 0, 0.985, e0)',
                overwrite='true',
                quiet='true')
        print 'Done!'
        
        print 'Calculating surface temperature (Ts) - K...'
        grass.mapcalc('Ts = $LS8_corre10/(1+((10.8*$LS8_corre10)/14380)*log($eNBf))',
              LS8_corre10='LS8_corre10', 
              eNBf='eNBf',               
              overwrite='true',
              quiet='true')
        print 'Done!'
        
        Ts_median=g.parse_command('r.univar', flags='ge', map='Ts', quiet = True)['median']        
        for line in open(MTLfile[0]):
                if 'EARTH_SUN_DISTANCE' in line:
                        d=float(line.split('=')[-1])
                elif 'RADIANCE_MAXIMUM_BAND_2' in line:
                        RADIANCE_MAXIMUM_BAND_2=float(line.split('=')[-1])
                elif 'RADIANCE_MAXIMUM_BAND_3' in line:
                        RADIANCE_MAXIMUM_BAND_3=float(line.split('=')[-1])
                elif 'RADIANCE_MAXIMUM_BAND_4' in line:
                        RADIANCE_MAXIMUM_BAND_4=float(line.split('=')[-1])
                elif 'RADIANCE_MAXIMUM_BAND_5' in line:
                        RADIANCE_MAXIMUM_BAND_5=float(line.split('=')[-1])
                elif 'RADIANCE_MAXIMUM_BAND_6' in line:
                        RADIANCE_MAXIMUM_BAND_6=float(line.split('=')[-1])
                elif 'RADIANCE_MAXIMUM_BAND_7' in line:
                        RADIANCE_MAXIMUM_BAND_7=float(line.split('=')[-1])
                elif 'REFLECTANCE_MAXIMUM_BAND_2' in line:
                        REFLECTANCE_MAXIMUM_BAND_2=float(line.split('=')[-1])
                elif 'REFLECTANCE_MAXIMUM_BAND_3' in line:
                        REFLECTANCE_MAXIMUM_BAND_3=float(line.split('=')[-1])
                elif 'REFLECTANCE_MAXIMUM_BAND_4' in line:
                        REFLECTANCE_MAXIMUM_BAND_4=float(line.split('=')[-1])
                elif 'REFLECTANCE_MAXIMUM_BAND_5' in line:
                        REFLECTANCE_MAXIMUM_BAND_5=float(line.split('=')[-1])
                elif 'REFLECTANCE_MAXIMUM_BAND_6' in line:
                        REFLECTANCE_MAXIMUM_BAND_6=float(line.split('=')[-1])
                elif 'REFLECTANCE_MAXIMUM_BAND_7' in line:
                        REFLECTANCE_MAXIMUM_BAND_7=float(line.split('=')[-1])
                elif 'SUN_ELEVATION' in line:
                        SUN_ELEVATION=float(line.split('=')[-1])
        ESUN_B2=(math.pi*d*d)*(RADIANCE_MAXIMUM_BAND_2/(REFLECTANCE_MAXIMUM_BAND_2*math.cos(math.radians(SUN_ELEVATION)))) 
        ESUN_B3=(math.pi*d*d)*(RADIANCE_MAXIMUM_BAND_3/(REFLECTANCE_MAXIMUM_BAND_3*math.cos(math.radians(SUN_ELEVATION)))) 
        ESUN_B4=(math.pi*d*d)*(RADIANCE_MAXIMUM_BAND_4/(REFLECTANCE_MAXIMUM_BAND_4*math.cos(math.radians(SUN_ELEVATION)))) 
        ESUN_B5=(math.pi*d*d)*(RADIANCE_MAXIMUM_BAND_5/(REFLECTANCE_MAXIMUM_BAND_5*math.cos(math.radians(SUN_ELEVATION)))) 
        ESUN_B6=(math.pi*d*d)*(RADIANCE_MAXIMUM_BAND_6/(REFLECTANCE_MAXIMUM_BAND_6*math.cos(math.radians(SUN_ELEVATION)))) 
        ESUN_B7=(math.pi*d*d)*(RADIANCE_MAXIMUM_BAND_7/(REFLECTANCE_MAXIMUM_BAND_7*math.cos(math.radians(SUN_ELEVATION)))) 
        ESUN=[ESUN_B2,ESUN_B3,ESUN_B4,ESUN_B5,ESUN_B6,ESUN_B7]
        print ESUN
        W = []
        for i in range(len(ESUN)):
                W += [ESUN[i]/sum(ESUN)]
        W2=W[0]
        W3=W[1]
        W4=W[2]
        W5=W[3]
        W6=W[4]
        W7=W[5]
        print W
        
        print 'Calculating albedo at top of atmosphere (aTOA)...'
        grass.mapcalc('aTOA=$LS8_corre2*$W2+$LS8_corre3*$W3+$LS8_corre4*$W4+$LS8_corre5*$W5+$LS8_corre6*$W6+$LS8_corre7*$W7',
                      LS8_corre2='LS8_corre2',W2=W2, 
                      LS8_corre3='LS8_corre3',W3=W3, 
                      LS8_corre4='LS8_corre4',W4=W4, 
                      LS8_corre5='LS8_corre5',W5=W5, 
                      LS8_corre6='LS8_corre6',W6=W6, 
                      LS8_corre7='LS8_corre7',W7=W7,         
                      overwrite='true',
                      quiet='true')
        print 'Done!'
        
        print 'Median surface temperature:', Ts_median,'K'
        
        print 'Calculating shortwave transmissivity of air (Tsw)...'
        grass.mapcalc('Tsw=0.75+0.00002*$MDT_Sebal',
                      MDT_Sebal='MDT_Sebal',
                      overwrite='true',
                      quiet='true')
        print 'Done!'
        
        print 'Calculating surface albedo (aS)...'
        grass.mapcalc('aS=($aTOA-0.03)/($Tsw^2)',
                      aTOA='aTOA',
                      Tsw='Tsw',
                      overwrite='true',
                      quiet='true')
        print 'Done!'
        
        print 'Calculating incoming shortwave radiation (Rsi) - W/m2...'
        grass.mapcalc('Rsi=1367*cos(90-$SUN_ELEVATION)*(1/$d)^2*$Tsw',
                      SUN_ELEVATION=SUN_ELEVATION,
                      Tsw='Tsw',
                      d=d,
                      overwrite='true',
                      quiet='true')
        print 'Done!'
        
        print 'Calculating outgoing longwave radiation (RLo) - W/m2...'
        grass.mapcalc('RLo=$e0f*5.67e-8*$Ts^4',
                      Ts='Ts',
                      e0f='e0f',
                      overwrite='true',
                      quiet='true')
        print 'Done!'

print 'Making the cold pixel mask...'
Ts_median=g.parse_command('r.univar', flags='ge', map='Ts', quiet = True)['median']        
grass.mapcalc('Pcold=if($NDVI>0.4 && $Ts<$Ts_median,$Ts,null())',
              NDVI='NDVI',
              aS='aS',
              Ts='Ts',
              Ts_median=Ts_median,
              overwrite='true',
              quiet='true')
print 'Choose the cold pixel coordinates in irrigation areas. Use the raster Pcold and CC_652 for help...'
xy_Pcold = str(input('Place coordinates (east,north): ')).strip('()')
print xy_Pcold
TsPcold_z=g.parse_command('r.what', map='Ts', coordinates=xy_Pcold)
z_TsPcold=float(dict.keys(TsPcold_z)[0].split('|')[3])
print 'Cold pixel temperature:',z_TsPcold,'K'
grass.write_command('v.in.ascii', input='-', output='Pcold', sep=',', stdin=xy_Pcold, overwrite='true', quiet='true')
print 'Calculating incoming longwave radiation (RLi) - W/m2...'
grass.mapcalc('RLi=0.85*((-log($Tsw))^0.09)*5.67e-8*$z_TsPcold^4',
              z_TsPcold=z_TsPcold,
              Tsw='Tsw',
              overwrite='true',
              quiet='true')
print 'Done!'
print 'Calculating net radiation flux (Rn) - W/m2...'
grass.mapcalc('Rn=(1-$aS)*$Rsi+$RLi-$RLo-(1-$e0f)*$RLi',
              aS='aS',
              RLi='RLi',
              RLo='RLo',
              Rsi='Rsi',
              e0f='e0f',
              overwrite='true',
              quiet='true')
print 'Done!'
print 'Calculating soil heat flux (G) - W/m2...'
grass.mapcalc('G_Rn=if(NDVI<0,0.5,(($Ts-273.15)/$aS)*(0.0038*$aS+0.0074*$aS^2)*(1-0.98*$NDVI^4))',
              aS='aS',
              Ts='Ts',
              #Rn='Rn',
              NDVI='NDVI',
              overwrite='true',
              quiet='true')
grass.mapcalc('G=$G_Rn*$Rn',
              G_Rn='G_Rn',
              Rn='Rn',
              overwrite='true',
              quiet='true')
print 'Done!'
print 'Making the hot pixel mask...'
grass.mapcalc('Phot=if($SAVI>0.18 && $SAVI<0.3, $Ts,null())',
              SAVI='SAVI',
              aS='aS',
              Ts='Ts',
              Ts_median=Ts_median,
              #Ts_max=Ts_max,
              overwrite='true',
              quiet='true')
print 'Choose the hot pixel coordinates in bare soil areas. Use the raster Phot and CC_652 for help...'
xy_Phot = str(input('Place coordinates (east,north): ')).strip('()')
print xy_Phot
print 'Done!'
TsPhot_z=g.parse_command('r.what', map='Ts', coordinates=xy_Phot)
z_TsPhot=float(dict.keys(TsPhot_z)[0].split('|')[3])
print 'Hot pixel temperature:',z_TsPhot,'K'
grass.write_command('v.in.ascii', input='-', output='Phot', sep=',', stdin=xy_Phot, overwrite='true', quiet='true')
print 'Calculating friction velocity (u*) for weather station - m/s...'
h=0.15
Zom=0.123*h
u_ast=0.41*u_2m/(math.log(2/Zom))
u_200m=u_ast*(math.log(200/Zom))/0.41
print 'Done!'
print 'Friction velocity', u_ast, 'm/s'
print 'Calculating the momentum roughness length map (Z0map) - m...'
grass.mapcalc('Z0map=exp(-5.809+5.62*$SAVI)',
              SAVI='SAVI',
              overwrite='true',
              quiet='true')
print 'Done!'
print 'Calculating the friction velocity map (u*map) - m/s...'
grass.mapcalc('u_astmap=0.41*$u_200m/log(200/$Z0map)',
              Z0map='Z0map',
              u_200m=u_200m,
              overwrite='true',
              quiet='true')
print 'Done!'
print 'Calculating aerodynamic resistance to heat transport map in terms of neutral stability (rah) - s/m...'
grass.mapcalc('rah=log(2/0.1)/($u_astmap*0.41)',
              u_astmap='u_astmap',
              overwrite='true',
              quiet='true')
print 'Done!'

GPhot_z=g.parse_command('r.what', map='G', coordinates=xy_Phot)
z_GPhot=float(dict.keys(GPhot_z)[0].split('|')[3])
rahPhot_z=g.parse_command('r.what', map='rah', coordinates=xy_Phot)
z_rahPhot=float(dict.keys(rahPhot_z)[0].split('|')[3])
RnPhot_z=g.parse_command('r.what', map='Rn', coordinates=xy_Phot)
z_RnPhot=float(dict.keys(RnPhot_z)[0].split('|')[3])
print 'Aerodynamic resistance correction, be patient...'

z_rahPhot_i=0
i=1
while (abs(z_rahPhot_i - z_rahPhot) > 0.00001):
        print 'Iteration number:',i 
        i = i+1
        z_rahPhot_i=z_rahPhot
        print 'Equation coefficents estimate >> dT = a.Ts + b'
        a=(z_RnPhot-z_GPhot)* z_rahPhot_i/((z_TsPhot-z_TsPcold)*1.25*1004)
        b=-a*z_TsPcold
        print 'a:',a, 'b:', b 
        print 'Calculating dT map - K...'
        grass.mapcalc('dT=$a*$Ts+$b',
                      Ts='Ts',
                      a=a,
                      b=b,
                      overwrite='true',
                      quiet='true')
        print 'Done!'
        print 'Calculating sensible heat flux (H) - W/m2...'
        grass.mapcalc('H=($dT/$rah)*1.25*1004',
                      dT='dT',
                      rah='rah',
                      overwrite='true',
                      quiet='true')
        print 'Done!'
        print 'Calculating the Monin-Obukhov length map (L) - m...'
        grass.mapcalc('L=-(1.25*1004*$Ts*$u_astmap^3)/(0.41*9.81*$H)',
                      Ts='Ts',
                      u_astmap='u_astmap',
                      H='H',
                      overwrite='true',
                      quiet='true')
        print 'Done!'
        print 'Calculating atmospheric stability correction (L200m, L2m, L01m)...'
        grass.mapcalc('L200m=if($L<0,2*log((1+(1-16*(200/$L))^0.25)/2)+log((1+(1-16*(200/$L))^0.5)/2)-2*atan((1-16*(200/$L))^0.25)+0.5*3.14159265,if($L>0,-5*(2/$L),0))',
                      L='L',
                      overwrite='true',
                      quiet='true')
        grass.mapcalc('L2m=if($L<0,2*log((1+(1-16*(2/$L))^0.5)/2),if($L>0,-5*(2/$L),0))',
                      L='L',
                      overwrite='true',
                      quiet='true')
        grass.mapcalc('L01m=if($L<0,2*log((1+(1-16*(0.1/$L))^0.5)/2),if($L>0,-5*(0.1/$L),0))',
                      L='L',
                      overwrite='true',
                      quiet='true')
        print 'Done!'
        print 'Calculating corrected friction velocity map (u*map) - m/s...'
        grass.mapcalc('u_astmap=0.41*$u_200m/(log(200/$Z0map)-$L200m)',
                      Z0map='Z0map',
                      u_200m=u_200m,
                      L200m='L200m',
                      overwrite='true',
                      quiet='true')
        print 'Done!'
        print 'Calculating corrected aerodynamic resistance to heat transport (rah) - s/m...'
        grass.mapcalc('rah=(log(2/0.1)-$L2m+$L01m)/($u_astmap*0.41)',
                      u_astmap='u_astmap',
                      L2m='L2m',
                      L01m='L01m',
                      overwrite='true',
                      quiet='true')
        print 'Done!'
        rahPhot_z=g.parse_command('r.what', map='rah', coordinates=xy_Phot)
        z_rahPhot=float(dict.keys(rahPhot_z)[0].split('|')[3])

HPhot_z=g.parse_command('r.what', map='H', coordinates=xy_Phot)
z_HPhot=float(dict.keys(HPhot_z)[0].split('|')[3])
dTPhot_z=g.parse_command('r.what', map='dT', coordinates=xy_Phot)
z_dTPhot=float(dict.keys(dTPhot_z)[0].split('|')[3])
print 'Hot pixel results, verify if Rhhot - Ghot = Hhot'
print 'Hhot:', z_HPhot,
print 'rahhot:', z_rahPhot,
print 'Ghot:', z_GPhot,
print 'Rnhot:', z_RnPhot, 'dT:',z_dTPhot
print 'Calculating the latent heat flux (LET) - W/m2...'
grass.mapcalc('LET=Rn-G-H',
              H='H',
              Rn='Rn',
              G='G',
              overwrite='true',
              quiet='true')
print 'Done!'
print 'Calculating instantaneous evapotranspiration (ETi) - mm/h...'
grass.mapcalc('ETi=if(3600*(LET/(2.45*10^6))<0,0,3600*(LET/(2.45*10^6)))',
              LET='LET',
              overwrite='true',
              quiet='true')
print 'Done!'
print 'Calculating reference evapotranspiration fraction (ETof)...'
grass.mapcalc('ETof=$ETi/$EToi',
              ETi='ETi',
              EToi=EToi,
              overwrite='true',
              quiet='true')
print 'Done!'
print 'Calculating daily evapotranspiration (ETday) - mm/dia...'
grass.mapcalc('ETday=$ETof*$ETo',
              ETof='ETof',
              ETo=ETo,
              overwrite='true',
              quiet='true')
print 'Done!'

