# Programa Bolha T.
import string
import math
def main():
        print "programa bolha t" 
#seta a principais variaveis
        phi1 = 1
        phi2 = 1
        dphi1 = 1
        dphi2 = 1
        contador = 1
        i = 0
        sair = "nao"
        limite = 0.000000001
##variaveis de antoine
        a1 = 8.07131
        a2 = 7.43155
        b1 = 1730.63
        b2 = 1554.679
        c1 = 233.426
        c2 = 240.337
##var do gamma - van laar 
        a = 1.6742
        b = 2.5452
##var do poynting 
        vl1 = 18.07 * 10 ** (-6)
        vl2 = 85.71 * 10 ** (-6)
##constantes
        r = 8.314*10 ** (-5)
        p = 1.01308
##var do phi
        tc1 = 647.3
        tc2 = 587.0
        pc1 = 221.2
        pc2 = 52.1
        vc1 = 57.1 * 10 ** (-6) 
        vc2 = 238.0 * 10 ** (-6)
        zc1 = 0.235
        zc2 = 0.254
        w1 = 0.344
        w2 = 0.281
        k12 = 1-((8*(vc1*vc2)**0.5)/((vc1**(1.0/3.0)+vc2**(1.0/3.0))**3)) 
        w = []
        dw = []

  # seta os arquivos de entrada e saida  
        #infilename = raw_input("entre com o nome do arquivo de entrada (entrada.txt): ")
        infilename = ""
        outfilename = raw_input("entre com o nome do arquivo de saida (saida.txt): ")
        if infilename == "":
                infilename = "entrada.txt"
        if outfilename == "":
                outfilename = "saida.txt"
                
               
        # abre o arquivo para leitura e para escrita respectivamente
        infile = open(infilename, 'w')
        outfile = open(outfilename, 'w')
                
        def psat(t, dt):
            t = t - 273.15
            dt = dt - 273.15
            p1sat = 10 ** (a1-b1/(c1+t))
            p2sat = 10 ** (a2-b2/(c2+t))
            dp1sat = 10 ** (a1-b1/(c1+dt))
            dp2sat = 10 ** (a2-b2/(c2+dt))
                
            p1satbar = p1sat * 0.001333
            p2satbar = p2sat * 0.001333
            dp1satbar = dp1sat * 0.001333
            dp2satbar = dp2sat * 0.001333
            return p1satbar, p2satbar, dp1satbar, dp2satbar
            
        def gamma(x1, x2):
            """calculo de van laar 1 e 2"""
            gamma1 = math.exp(a / (1 + a * x1/(b * x2)) ** 2.0)            
            gamma2 = math.exp(b / (1 + b * x2/(a * x1)) ** 2.0)		
            return gamma1, gamma2
        

        def poynting(t, dt):
            """calculo da correcao de poynting"""
            p1satbar, p2satbar, dp1satbar, dp2satbar = psat(t, dt)
            poy1 = math.exp(vl1 * (p - p1satbar)/(r * t))
            poy2 = math.exp(vl2 * (p - p2satbar)/(r * t))	
            dpoy1 = math.exp(vl1 * (p - dp1satbar)/(r * dt))
            dpoy2 = math.exp(vl2 * (p - dp2satbar)/(r * dt))
            return poy1, poy2, dpoy1, dpoy2
            
        def phisat(t, dt):
            """calculo de phisat"""
            tr1 = t/tc1
            tr2 = t/tc2
            dtr1 = dt/tc1
            dtr2 = dt/tc2
        
            b10= 0.083-0.422/(tr1**1.6)
            b20= 0.083-0.422/(tr2**1.6)
            db10= 0.083-0.422/(dtr1**1.6)
            db20= 0.083-0.422/(dtr2**1.6)
            b11 = 0.139-0.172/(tr1**4.2)
            b21 = 0.139-0.172/(tr2**4.2)
            db11 = 0.139-0.172/(dtr1**4.2)
            db21 = 0.139-0.172/(dtr2**4.2)
            B11 = r*tc1/pc1 * (b10+w1*b11)
            B22 = r*tc2/pc2 * (b20+w2*b21)
            dB11 = r*tc1/pc1 * (db10+w1*db11)
            dB22 = r*tc2/pc2 * (db20+w2*db21)
            
            p1satbar, p2satbar, dp1satbar, dp2satbar = psat(t, dt)
            phi1sat = math.exp(B11 * p1satbar/(r*t))
            phi2sat = math.exp(B22 * p2satbar/(r*t))
            dphi1sat = math.exp(dB11 * dp1satbar/(r*dt))
            dphi2sat = math.exp(dB22 * dp2satbar/(r*dt))
            return phi1sat, phi2sat, dphi1sat, dphi2sat 
            
        def calculoy(t, dt):
            """calculo do y"""
            p1satbar, p2satbar, dp1satbar, dp2satbar = psat(t, dt)
            y1 = x1 * gamma1 * p1satbar * phi1sat * poy1/ (phi1 * p)
            y2 = x2 * gamma2 * p2satbar * phi2sat * poy2 / (phi2 * p)
            dy1 = x1 * gamma1 * dp1satbar * dphi1sat * dpoy1/ (dphi1 * p)
            dy2 = x2 * gamma2 * dp2satbar * dphi2sat * dpoy2 / (dphi2 * p)
            return y1, y2, dy1, dy2
            
        def phi(t, dt, y1, y2, dy1, dy2):
            """calculo do phi"""
            tr1 = t/tc1
            tr2 = t/tc2
            dtr1 = dt/tc1
            dtr2 = dt/tc2
            tc12 = (tc1*tc2) ** 0.5 * (1-k12)
            vc12 = ((vc1 ** (1.0/3.0) + vc2 ** (1.0/3.0))/2.0)**3.0
            zc12 = (zc1+zc2)/2.0
            pc12 = zc12*r*tc12/vc12
            w12 = (w1+w2)/2.0
            tr12 = t/tc12
            dtr12 = dt/tc12
        
            b10= 0.083-0.422/(tr1**1.6)
            b20= 0.083-0.422/(tr2**1.6)
            b120= 0.083-0.422/(tr12**1.6)
            db10= 0.083-0.422/(dtr1**1.6)
            db20= 0.083-0.422/(dtr2**1.6)
            db120= 0.083-0.422/(dtr12**1.6)

            b11 = 0.139-0.172/(tr1**4.2)
            b21 = 0.139-0.172/(tr2**4.2)
            b121 = 0.139-0.172/(tr12**4.2)
            db11 = 0.139-0.172/(dtr1**4.2)
            db21 = 0.139-0.172/(dtr2**4.2)
            db121 = 0.139-0.172/(dtr12**4.2)

            B11 = r*tc1/pc1 * (b10+w1*b11)
            B22 = r*tc2/pc2 * (b20+w2*b21)
            B12 = r*tc12/pc12 * (b120+w12*b121)
            dB11 = r*tc1/pc1 * (db10+w1*db11)
            dB22 = r*tc2/pc2 * (db20+w2*db21)
            dB12 = r*tc12/pc12 * (db120+w12*db121)

            b = y1**2*B11+2*y1*y2*B12+y2**2*B22
            db = dy1**2*dB11+2*dy1*dy2*dB12+dy2**2*dB22
            #~ print (x2)
            #~ print (y1)
            #~ print (B11)
            #~ print (y2)
            #~ print (B12) 
            #~ print (b) 
            #~ print (p) 
            #~ print (r) 
            #~ print (t)           
            phi1 = math.exp((2.0*y1*B11+2.0*y2*B12-b)*p/(r*t))
            phi2 = math.exp((2.0*y2*B22+2.0*y1*B12-b)*p/(r*t))
            dphi1 = math.exp((2.0*dy1*dB11+2.0*dy2*dB12-db)*p/(r*dt))
            dphi2 = math.exp((2.0*dy2*dB22+2.0*dy1*dB12-db)*p/(r*dt))
            return phi1, phi2, dphi1, dphi2
            
        def ajustaT(t, dt):
            """Ajuste de T"""
            p1satbar, p2satbar, dp1satbar, dp2satbar = psat(t, dt)            
            #~ gamma1, gamma2 = gamma(x1, x2)
            #~ poy1, poy2, dpoy1, dpoy2 = poynting(t, dt)
            #~ phi1sat, phi2sat, dphi1sat, dphi2sat = phisat(t, dt)
            #~ y1, y2, dy1, dy2 = calculoy(t, dt)
            #~ phi1, phi2, dphi1, dphi2 = phi(t, dt, y1, y2, dy1, dy2)                
            
            ft = 1-y1_ant*phi1*p/(gamma1*p1satbar*phi1sat*poy1)-y2_ant*phi2*p/(gamma2*p2satbar*phi2sat*poy2)
            fdt = 1-dy1_ant*dphi1*p/(gamma1*dp1satbar*dphi1sat*dpoy1)-dy2_ant*dphi2*p/(gamma2*dp2satbar*dphi2sat*dpoy2)
            dft = (ft-fdt)/0.005
            tnova = t-ft/dft
            t = tnova
            return t
        
        # Processa cada linha do arquivo de entrada
        while str(sair) <> "sim":            
            #temperatura, x1 = string.split(line, ";")
            # pega a o x1 
            x1 = raw_input("Entre com o valor de x1: ")
            #grava no arquivo entrada.txt o x1
            infile.write(str(x1) +'\n')
            #transforma em float
            x1 = float(x1)     
            # recebe o valor de temperatura kelvin
            t = raw_input("estime o valor de temperatura em kelvin: ")
            t = float(t)
            dt = t - 0.005
            x2 = float(1 - x1)  
            # Chamada de Funcoes            
            gamma1, gamma2 = gamma(x1, x2)
            poy1, poy2, dpoy1, dpoy2 = poynting(t, dt)
            phi1sat, phi2sat, dphi1sat, dphi2sat = phisat(t, dt)
            y1, y2, dy1, dy2 = calculoy(t, dt)
            w.append (y1 + y2)
            dw.append(dy1 + dy2)            
            
            if contador == 1:    
                y1 = y1 / w[i]
                y2 = y2 / w[i]                
                dy1 = dy1 / dw[i]
                dy2 = dy2 / dw[i]
                #phi1 e phi 2 e o contador para o proximo x deve ser 1
                phi1, phi2, dphi1, dphi2 = phi(t, dt, y1, y2, dy1, dy2)                                  
                y1, y2, dy1, dy2 = calculoy(t, dt)                
                w.append(y1 + y2)
                dw.append(dy1 + dy2)
                contador = contador + 1 
            #~ print w[i]
            #~ print w[i+1]
            while round(w[i], 8)!= 1:                
                while (abs(w[i+1] - w[i]) >= limite):
                    #normaliza                    
                    y1 = y1 / w[i+1]
                    #guarda o valor anterior
                    y1_ant = y1                   
                    y2 = y2 / w[i+1]
                    #guarda o valor anterior
                    y2_ant = y2                    
                    dy1 = dy1 / dw[i+1]
                    #guarda o valor anterior
                    dy1_ant = dy1                    
                    dy2 = dy2 / dw[i+1]
                    #guarda o valor anterior
                    dy2_ant = dy2
                    #
                    phi1, phi2, dphi1, dphi2 = phi(t, dt, y1, y2, dy1, dy2)
                    y1, y2, dy1, dy2 = calculoy(t, dt)                   
                    w.append(y1 + y2)
                    dw.append(dy1 + dy2)
                    i = i + 1                  
                    #~ dw[i+2] =  dy1 + dy2
                    #~ w[i+2] = y1 + y2
                # o dt muda sempre e no incio tem q ter valor
                t = ajustaT(t, dt)
                dt = t - 0.005
                gamma1, gamma2 = gamma(x1, x2)
                POY1, POY2, dPOY1, dPOY2 = poynting(t, dt)
                phi1sat, phi2sat, dphi1sat, dphi2sat = phisat(t, dt)
                y1, y2, dy1, dy2 = calculoy(t, dt)
                w.append(y1 + y2)
                dw.append(dy1 + dy2)
                i = i + 1 
                #~ print w[i]
                #~ print w[i+1]
                print "t: " + str(t) +'\n'
                print "y1: " +str(y1) +'\n'
                #~ w[i+1] =  y1 + y2     
                #~ dw[i+1] = dy1 + dy2
            
            # Escreve no arquivo de saida o x2
            outfile.write(str(x2) +'\n')
            #Retorna o valor das var para o proximo x1
            phi1 = 1
            phi2 = 1
            dphi1 = 1
            dphi2 = 1
            contador = 1
            i = 0
            #Verifica se deseja sair
            sair = raw_input("Deseja sair? 'sim', 'nao': ")
            

# close both files
        infile.close()
        outfile.close()
        print "Arquivo de saida com o valor de X2 criado: ", outfilename
main()
