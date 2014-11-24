#!/usr/bin/env ruby

def gen_plot_cmd(nbProc, nbParticle)
	color = ["red", "green", "blue", "orange", "purple"]
	cmd = "plot "
	for i in 1..nbProc
		s= 2
		e = 3
		for j in 1..nbParticle
			cmd += "\'datafile#{i-1}\' index (i-1) using #{s}:#{s+1} with circle lc rgb \'#{color[i-1% color.size]}\' title 'Processeur#{i-1}, Particule#{j-1}'"
			s+=2
			unless (i== nbProc && j==nbParticle)
				cmd+= ', '
			end
		end
	end
	cmd
end

def print_on_file(filename, s)
	system("echo \"#{s}\" " + ">> #{filename}")	
end

def gen_plot_script(filename, output, nbProc, nbParticle)
	File.open(filename, 'w+') do |f|  
		f.puts("set terminal gif size 800,600 animate delay 10\n"+
				"set output '#{output}'\n"+
				"stats 'datafile0' nooutput\n"+
				"set xrange [-2000:2001]\n"+
				"set yrange [-2000:2001]\n"+
				"do for [i=1:int(STATS_blocks)] {\n\t"+
				gen_plot_cmd(nbProc, nbParticle)+
				"\n}\n")
	end 
end


def simulate(nbProcess, nbParticle, pattern)
	filename = "particles#{nbProcess}#{nbParticle}#{pattern}"
	plotfile = filename+".plt"
	giffile = filename+".gif"
	system("mpiexec -np #{nbProcess} ./simulation.out #{pattern}")
	gen_plot_script(plotfile, giffile, nbProcess, nbParticle)
	system("gnuplot #{plotfile}")
	system("eog #{giffile}")
end

#simulate(2, 2, "random")

gen_plot_script("particles", 'a.plt', 4,1)