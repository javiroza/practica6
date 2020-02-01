# Format i nom de la imatge
set term png
set output "P6-1920-fig3.png"

# Mostra els eixos
#set xzeroaxis
#set yzeroaxis

# Títol del gràfic
set title "Convergència de I_3"

# Rang dels eixos
#set xrange[-1.00e-10:1.00e10]
#set yrange[-1.00e-10:1.00e10]

# Títols dels eixos
set xlabel "Nombre d'iteracions, N (U.A.)"
set ylabel "Valor de la Integral (U.A.)"

# Canvia els nombres dels eixos per nombres personalitzats
#set ytics("1x10^-^1^0" 1.00e-10,"1x10^-^0^5" 1.00e-05,"1x10^0" 1.00e+00,"1x10^5" 1.00e+05,"1x10^1^0" 1.00e+10)
#set xtics("1x10^-^3" 1.00e-03,"1x10^-^2" 1.00e-02,"1x10^-^1" 1.00e-01,"1x10^0" 1.00e+00,"1x10^1" 1.00e+01)

# Format dels nombres dels eixos
set format y '%.4f'
#set format x '%.1e'

# Escala dels eixos logarítmica
#set logscale y
#set logscale x

# Posició de la llegenda
set key top right

# Plot 
plot "P6-1920-res.dat" index 2 using 1:2 with points t "I_3(N)"
#pause -1