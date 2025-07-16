
# Exemplo de execucao:
# julia trigger_arc_tsp.jl --inputfilename ../instances/instances_release_1/grf1.txt \
#  --ra 123456 --maxtime_ilp 100 --maxtime_lb_lp 100 \
#  --maxtime_ub_lp 100
  
#--------------------------------------------------------
#        Nao mude o conteudo deste arquivo
#--------------------------------------------------------
using ArgParse
using Base
using Random

#--------------------------------------------------------
# Mude apenas o conteudo do arquivo trigger_arc_tsp_routines.jl

include("trigger_arc_tsp_include.jl")	# Nao mude este arquivo
include("trigger_arc_tsp_routines.jl")	# Coloque suas rotinas neste arquivo
include("exemplo_usando_pli.jl")
include("helper.jl")

function main()
	 
	(inputfilename,	# input file as given in the internet
	seednumber,		# Randomized steps become deterministic for fixed seed
	maxtime_lb_lp,		# Max. time to compute LB via LP (may use cutting planes)
	maxtime_lb_rlxlag,	# Max. time to compute LB by Lagrangean Relaxation
	maxtime_lb_colgen,	# Max. time to compute LB via LP column generation
	maxtime_ub_lp,	# Max. time to compute UB via LP rounding
	maxtime_ub_rlxlag,	# Max. time to compute UB via Lagrangean Relaxation
	maxtime_ub_colgen,  	# Max. time to compute UB via LP column generation
	maxtime_ilp,		# Max. time to compute LB and UB via exact Branch and Cut
	ra,
	logfilename) = getparameters()
	
	println()
	println("---------- PARAMETER VALUES ... -----------")
	println("inputfilename =      ",inputfilename)
		println("seednumber =         ",seednumber)
	println("maxtime_lb_lp =      ",maxtime_lb_lp)
	println("maxtime_lb_rlxlag =  ",maxtime_lb_rlxlag)
	println("maxtime_lb_colgen =  ",maxtime_lb_colgen)
	println("maxtime_ub_lp =      ",maxtime_ub_lp)
	println("maxtime_ub_rlxlag =  ",maxtime_ub_rlxlag)
	println("maxtime_ub_colgen =  ",maxtime_ub_colgen)
	println("maxtime_ilp =        ",maxtime_ilp,)
	println("ra =                 ",ra)

	T=TriggerArcTSP(inputfilename,seednumber,
			maxtime_lb_lp,
			maxtime_lb_rlxlag, 
			maxtime_lb_colgen,
			maxtime_ub_lp,
			maxtime_ub_rlxlag,
			maxtime_ub_colgen,
			maxtime_ilp,
			ra,
			logfilename)

	if(T.NNodes == 0)
		println("Instance is too big. Skiping...")
		return
	end

	 # # O proximo trecho permite verificar se o arquivo lido esta armazenado
	 # # corretamente. O trecho escreve os dados para um novo arquivo e
	 # # compara com diff.
	 # # Write the instance again to a text file, so that we can compare with diff
	 # # The next lines are valid only if the corresponding parts of the inputfile
	 # # corresponding to arcs and triggers are given in the order of arc and 
	 # # trigger id's. (there are some instances that do not follow the id's order)
	 # outputfilename=inputfilename*"_copy"
	 # write_instance(T,outputfilename)
	 # if (isfile(inputfilename)) &&  (isfile(outputfilename)) 
	 #    output=read(`diff --brief -s $inputfilename $outputfilename`, String)
   	 #    println(output)
	 # end

	# print_instance(T)

	#TriggerArcTSP_lb_lp(T)
	#WriteLogFile(T,"lb_lp")
	
	# TriggerArcTSP_lb_rlxlag(T)
	# WriteLogFile(T,"lb_rlxlag")
	
	# TriggerArcTSP_lb_colgen(T)
	#WriteLogFile(T,"lb_colgen")
	
	# TriggerArcTSP_ub_lp(T)
	# WriteLogFile(T,"ub_lp")

	# TriggerArcTSP_ub_rlxlag(T)
	# WriteLogFile(T,"ub_rlxlag")

	# TriggerArcTSP_ub_colgen(T)
	#WriteLogFile(T,"ub_colgen")
	
	TriggerArcTSP_ilp(T)
	WriteLogFile(T,"ilp")

	# Exemplo_PLI(T,10)
end


# --------------------------------------------------------------------
main()
