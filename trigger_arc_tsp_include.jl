# --------------------------------------------------------
#        Nao mude o conteudo deste arquivo
# --------------------------------------------------------


# ------------------------------------------------------------
# Exemplo de execucao:
#
# julia tsp.jl --time 5 --nnodes 10 --seed 123 --distfactor 0.5
#
# time: tempo maximo que seu programa pode executar
# nnodes: numero de no's do grafo
# seed: semente para iniciar a geracao de numeros aleatorios
# 	num mesmo programa, para testar com instancias diferentes, basta
#	usar outro numero aleatorio.
# distfactor: para gerar o numero aleatorio, a entrada usa a distancia euclidiana
# 	e acrescenta mais um fator percentual aleatorio dentro deste fator.
#
# --------------------------------------------------------
#        Nao mude o conteudo deste arquivo
# --------------------------------------------------------
# Se for primeira vez, precisa instalar pacotes:
# using Pkg
# Pkg.add("ArgParse")
# Pkg.add("Random")
# Pkg.add("Parameters")

using ArgParse
using Random
using Parameters

#--------------------------------------------------------

const IntType   = Int64
const MaxInt	= typemax(IntType)
const MinInt	= typemin(IntType)

const FloatType = Float64
const MaxFloat	= typemax(FloatType)
const MinFloat	= typemin(FloatType)

# ===========================================
mutable struct NodeType
    NodeName ::String
    px ::IntType # if we know the node position, we can draw the graph/solution
    py ::IntType
    # -------------------------------------------
    function NodeType()
    	     NodeType("",-1,-1)
    	     end
    function NodeType(n::IntType)
    	     NodeType(string(n),-1,-1)
    	     end
    function NodeType(n::IntType, px::IntType, py::IntType)
    	     NodeType(string(n),px,py)
    	     end
    function NodeType(a::String, px::IntType, py::IntType)
    	     new(a,px,py)
	     end
end

# ===========================================
mutable struct ArcType
    	u ::IntType  # Index of 'from' node
    	v ::IntType  # Index of 'to' node
    	cost ::FloatType # Cost of the arc

    	function ArcType()
    	     	 ArcType(-1,-1,0.0)
    	     	 end
    	function ArcType(uu::IntType, vv::IntType, cc::FloatType)
    	     	 new(uu,vv,cc)
	     	 end
end
# ===========================================
mutable struct TriggerType
	trigger_arc_id::IntType
	target_arc_id::IntType
	cost::FloatType # new target cost, in case activated by trigger_arc
    	function TriggerType()
    	     	 TriggerType(-1,-1,0.0)
    	     	 end
    	function TriggerType(a_trigger::IntType, a_target::IntType, cc::FloatType)
    	     	 new(a_trigger,a_target,cc)
	     	 end
end

# # ===========================================
# mutable struct SolutionType
# 	SolName::String         # Name of the method used to obtain the solution
# 	SolArc::Vector{IntType} # Solution is a sequence of arc indexes.
# 				# No. of arcs in SolArc is given by the vector size
# 	function SolutionType()
# 		 SolutionType("", MaxFloat, IntType[])
# 	end
# 	function SolutionType(  _SolName::String,
# 				_SolValue::FloatType,
# 				_SolArc::Vector{IntType})
# 		 new(_SolNArcs,_SolValue,_SolArc)
# 	end
# end

# ===========================================
mutable struct TriggerArcTSP
	# Info about the input
	NNodes::IntType
	NArcs::IntType
	NTriggers::IntType
	Node::Vector{NodeType}
	Arc::Vector{ArcType}
	Trigger::Vector{TriggerType}

	# Command line parameters 
	inputfilename::String	# input file as given in the internet
	seednumber::IntType	# Randomized routines become deterministic for fixed seed

	maxtime_lb_lp::IntType       # Max. time to compute LB by LP (may use cut.planes)
	maxtime_lb_rlxlag::IntType   # Max. time to compute LB by Lagran. Relax.
	maxtime_lb_colgen::IntType   # Max. time to compute LB via LP column gener.
	
	maxtime_ub_lp::IntType       # Max. time to compute UB by LP rounding
	maxtime_ub_rlxlag::IntType   # Max. time to compute UB by Lagran. Relax.
	maxtime_ub_colgen::IntType   # Max. time to compute UB via LP column gener.

	maxtime_ilp::IntType	     # Max. time to compute LB/UB by exact Branch and Cut

	# Output info: time used in each process
	time_lb_lp::IntType       # Time used to compute LB by LP (may use cut.planes)
	time_lb_rlxlag::IntType   # Time used to compute LB by Lagran. Relax.
	time_lb_colgen::IntType   # Time used to compute LB via LP column gener.

	time_ub_lp::IntType  	  # Time used to compute UB by LP rounding
	time_ub_rlxlag::IntType   # Time used to compute UB by Lagran. Relax.
	time_ub_colgen::IntType   # Time used to compute UB via LP column gener.

	time_ilp::IntType   	  # time used to compute LB/UB by exact Branch and Cut

	# Output info: computed lower bounds
	lb_lp::FloatType     # LB by LP (may use cutt. planes)
	lb_rlxlag::FloatType # LB by Lagrangean Relaxation
	lb_colgen::FloatType # LB via LP column gener. (OPTIONAL)
	lb_ilp::FloatType    # Best LB of exact Branch and Cut (equals ub_ilp if optimum)

	# Output info: value of the UB solutions obtained by each LP approach
	ub_lp::FloatType      # UB by LP rounding (value of heuristic/metaheuristic sol.)
	ub_rlxlag::FloatType  # UB by Lag. Relax. (value of heuristic/metaheuristic sol.)
	ub_colgen::FloatType  # UB via LP column gener. (heur. or metaheuristic based)
	ub_ilp::FloatType     # UB value of the best solution of exact Branch and Cut
	
	# Output info: set of arcs in the UB solutions obtained by each LP approach
	ub_lp_arcs::Vector{IntType}   	 # UB sol. by LP rounding (set of arcs)
	ub_rlxlag_arcs::Vector{IntType}  # UB sol. by Lag. Relax. (set of arcs)
	ub_colgen_arcs::Vector{IntType}  # UB sol. by column gener. (set of arcs)
	ub_ilp_arcs::Vector{IntType}  	 # best UB sol. of exact Branch and Cut (set of arcs)

	# Output info: other info interesting for the exact ILP algorithm
	nn_ilp::IntType		# Number of nodes used in the branch tree algorithm

	# User information
	ra::IntType

	# ============================================================================
	# Ordem dos campos:
	# NNodes,	  NArcs,	      NTriggers,
	# Node,		  Arc,		      Trigger, 
	# inputfilename,  seednumber,
	# maxtime_lb_lp,  maxtime_lb_rlxlag,  maxtime_lb_colgen,
	# maxtime_ub_lp,  maxtime_ub_rlxlag,  maxtime_ub_colgen,  maxtime_ilp 
	# time_lb_lp, 	  time_lb_rlxlag,     time_lb_colgen,
	# time_ub_lp,     time_ub_rlxlag,     time_ub_colgen,	  time_ilp
	# lb_lp,      	  lb_rlxlag, 	      lb_colgen, 	  lb_ilp,
	# ub_lp, 	  ub_rlxlag, 	      ub_colgen, 	  ub_ilp
	# ub_lp_arcs,     ub_rlxlag_arcs,     ub_colgen_arcs,	  ub_ilp_arcs
	# nn_ilp,	  ra
		 
	# ============================================================================
	# Constructors
	function TriggerArcTSP() # Constructor with no previous information
		 # Parameters must follow the field order 
		 TriggerArcTSP(
		 0,0,0,					# default for NNodes, NArcs, NTriggers
		 NodeType[],ArcType[],TriggerType[],	# default are empty vectors
		 "",1,		# default for inputfilename, seednumber
		 
		 MaxInt, 	# default for maxtime_lb_lp
		 MaxInt, 	# default for maxtime_lb_rlxlag
		 MaxInt, 	# default for maxtime_lb_colgen
		 
		 MaxInt, 	# default for maxtime_ub_lp
		 MaxInt, 	# default for maxtime_ub_rlxlag
		 MaxInt, 	# default for maxtime_ub_colgen
		 MaxInt, 	# default for maxtime_ilp to compute UB and LB
		 
		 0, 		# default for time_lb_lp
		 0,		# default for time_lb_rlxlag
		 0, 		# default for time_lb_colgen
		 
		 0, 		# default for time_ub_lp
		 0, 		# default for time_ub_rlxlag
		 0, 		# default for time_ub_colgen
		 0, 		# default for time_ilp

		 MinFloat,	# default for lb_lp
		 MinFloat,	# default for lb_rlxlag
		 MinFloat,	# default for lb_colgen
		 MinFloat,	# default for lb_ilp
		 
		 MaxFloat,	# default for ub_lp
		 MaxFloat,	# default for ub_rlxlag
		 MaxFloat,	# default for ub_colgen
		 MaxFloat,	# default for ub_ilp

		 Vector{IntType}(), # set of arcs in the ub_lp solution starts by an empty vector.
		 Vector{IntType}(), # set of arcs in the ub_rlxlag  solution starts by an empty vector.
		 Vector{IntType}(), # set of arcs in the ub_colgen  solution starts by an empty vector.
		 Vector{IntType}(), # set of arcs in the ub_ilp     solution starts by an empty vector.
		 
		 0,		# default for nn_ilp, the number of branch and bound/cut tree
		 999999)	# academic number
	end
		 
	function TriggerArcTSP( _inputfilename::String,
				_seednumber::IntType,
				_maxtime_lb_lp::IntType, 
			 	_maxtime_lb_rlxlag::IntType,
	 		 	_maxtime_lb_colgen::IntType,
				
			 	_maxtime_ub_lp::IntType,
				_maxtime_ub_colgen::IntType,
				_maxtime_ub_rlxlag::IntType,
			 	_maxtime_ilp::IntType,
				_ra::IntType)
		 fp=open(_inputfilename,"r")
		 header = readline(fp)
		 sNNodes,sNArcs,sNTriggers = split(header," ")
		 
		 _NNodes = parse(IntType,sNNodes)
		 _NArcs = parse(IntType,sNArcs)
		 _NTriggers = parse(IntType,sNTriggers)
		 
		 _Node::Vector{NodeType} = Vector{NodeType}(undef,_NNodes)
		 _Arc::Vector{ArcType} = Vector{ArcType}(undef,_NArcs)
		 _Trigger::Vector{TriggerType} = Vector{TriggerType}(undef,_NTriggers)
		 
		 for n in 1:_NNodes
		     _Node[n] = NodeType(string(n), rand(1:10000), rand(1:10000))
		 end
		 for a in 1:_NArcs
		     _Arc[a] = ArcType(-1, -1, 0.0)
		 end
		 for t in 1:_NTriggers
		     _Trigger[t] = TriggerType( -1, -1, 0.0 )
		 end
		 for a in 1:_NArcs
		     arcline = split(readline(fp)," ")
		     # println(arcline)
		     arc_id = parse(IntType,arcline[1])+1
		     arc_u = parse(IntType,arcline[2])+1
		     arc_v = parse(IntType,arcline[3])+1
		     arc_cost = parse(FloatType,arcline[4])
		     if (_Arc[arc_id].u != -1)
		        error("The arc ",a," is duplicated.")
		     else
		     	_Arc[arc_id].u = arc_u
		     end
		     if (_Arc[arc_id].v != -1)
		        error("The arc ",a," is duplicated.")
		     else
		     	_Arc[arc_id].v = arc_v
		     end
		     _Arc[arc_id].cost = arc_cost
		 end
		 # Verify if all arcs were read
 		 for a in 1:_NArcs
		     if (_Arc[a].u == -1) || (_Arc[a].v == -1)
		     	error("The arc ",a," is not defined.")
		     end
		 end
		 for t in 1:_NTriggers
		     tline = split(readline(fp)," ")
		     t_id = parse(IntType,tline[1])+1
		     if (t_id<1) || (t_id>_NTriggers)
		     	error("Trigger id outside valid range.")
		     end
		     # trigger arc
		     t_arctrig_id = parse(IntType,tline[2])+1
		     t_arctrig_u = parse(IntType,tline[3])+1
		     t_arctrig_v = parse(IntType,tline[4])+1
		     # target arc		     
		     t_arctarg_id = parse(IntType,tline[5])+1
		     t_arctarg_u = parse(IntType,tline[6])+1
		     t_arctarg_v = parse(IntType,tline[7])+1
		     t_cost = parse(FloatType,tline[8])

		     if ((_Trigger[t_id].trigger_arc_id != -1) ||
		     	(_Trigger[t_id].target_arc_id != -1))
		     	error("Trigger id appears duplicated.")
		     end
		     _Trigger[t_id].trigger_arc_id = t_arctrig_id
		     _Trigger[t_id].target_arc_id = t_arctarg_id
		     _Trigger[t_id].cost = t_cost
		     if ((_Arc[_Trigger[t_id].trigger_arc_id].u != t_arctrig_u) ||
		     	(_Arc[_Trigger[t_id].trigger_arc_id].v != t_arctrig_v) ||
		     	(_Arc[_Trigger[t_id].target_arc_id].u != t_arctarg_u) ||
		     	(_Arc[_Trigger[t_id].target_arc_id].v != t_arctarg_v))
			error("Trigger relation has wrong values.")
		     end
		 end
		 close(fp)

		 # We can represent a solution, say ub, as a vector indexed by the node ids.
		 # Given a node i (i in 1,...,NNodes) and an arc a=(i,j),
		 # we can have ub[i]=a. In this case, the arc living node i in the solution ub
		 # uses arc a. So, the next node after i is the node Arc[a].v
		 _ub_lp_arcs	 = Vector{IntType}(undef,_NNodes)
		 _ub_rlxlag_arcs = Vector{IntType}(undef,_NNodes)
		 _ub_colgen_arcs = Vector{IntType}(undef,_NNodes)
		 _ub_ilp_arcs    = Vector{IntType}(undef,_NNodes)
		 for i in 1:_NNodes
		     _ub_lp_arcs[i]     = -1
		     _ub_rlxlag_arcs[i] = -1
		     _ub_colgen_arcs[i] = -1
		     _ub_ilp_arcs[i]    = -1
		 end
		 
		 new(	_NNodes,_NArcs,_NTriggers,
			_Node,_Arc,_Trigger,
			_inputfilename,_seednumber,
			_maxtime_lb_lp,_maxtime_lb_rlxlag,_maxtime_lb_colgen,
			_maxtime_ub_lp,_maxtime_ub_rlxlag,_maxtime_ub_colgen,_maxtime_ilp,
			
		 	0, 		# default for time_lb_lp
		 	0,		# default for time_lb_rlxlag
		 	0, 		# default for time_lb_colgen
			
		 	0, 		# default for time_ub_lp
		 	0, 		# default for time_ub_rlxlag
		 	0, 		# default for time_ub_colgen
		 	0, 		# default for time_ilp
		 
		 	MinFloat,	# default for lb_lp
		 	MinFloat,	# default for lb_rlxlag
		 	MinFloat,	# default for lb_colgen
		 	MinFloat,	# default for lb_ilp
			
		 	MaxFloat,	# default for ub_lp
		 	MaxFloat,	# default for ub_rlxlag
		 	MaxFloat,	# default for ub_colgen
		 	MaxFloat,	# default for ub_ilp
			
		 	_ub_lp_arcs,
			_ub_rlxlag_arcs,
			_ub_colgen_arcs,
			_ub_ilp_arcs,

			0,		# default for nn_ilp
			_ra)		# academic number
	end
end

function write_instance(T::TriggerArcTSP, outputfilename::String)
	 fp=open(outputfilename,"w")
	 
	 # Write the header
	 println(fp,T.NNodes," ",T.NArcs," ",T.NTriggers)
	 
	 # Write the arcs
	 for a in 1:T.NArcs
	     println(fp,a-1," ",T.Arc[a].u-1," ",T.Arc[a].v-1," ",T.Arc[a].cost)
	 end
	 
	 # Write the Triggers
	 for t in 1:T.NTriggers
	     println(fp,t-1," ",
	     T.Trigger[t].trigger_arc_id-1," ",
	     T.Arc[T.Trigger[t].trigger_arc_id].u-1," ",
	     T.Arc[T.Trigger[t].trigger_arc_id].v-1," ",
	     T.Trigger[t].target_arc_id-1," ",
	     T.Arc[T.Trigger[t].target_arc_id].u-1," ",
	     T.Arc[T.Trigger[t].target_arc_id].v-1," ",
	     T.Trigger[t].cost)
	 end
	 close(fp)
end

function print_instance(T::TriggerArcTSP)
	 
	 # Print the header
	 println(T.NNodes," ",T.NArcs," ",T.NTriggers)
	 
	 # Print the arcs
	 for a in 1:T.NArcs
	     println(a-1," ",T.Arc[a].u-1," ",T.Arc[a].v-1," ",T.Arc[a].cost)
	 end
	 
	 # Print the Triggers
	 for t in 1:T.NTriggers
	     println(t-1," ",
	     T.Trigger[t].trigger_arc_id-1," ",
	     T.Arc[T.Trigger[t].trigger_arc_id].u-1," ",
	     T.Arc[T.Trigger[t].trigger_arc_id].v-1," ",
	     T.Trigger[t].target_arc_id-1," ",
	     T.Arc[T.Trigger[t].target_arc_id].u-1," ",
	     T.Arc[T.Trigger[t].target_arc_id].v-1," ",
	     T.Trigger[t].cost)
	 end
end



function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--inputfilename"
            help = "filename of the input instance"
	    arg_type = String
	    required = true
        "--seednumber"
            help = "seed for the random number generator"
	    arg_type = Int
	    required = false
	    default = 1
        "--maxtime_lb_lp"
            help = "maximum time, in seconds, to compute LB by LP Relaxation"
            arg_type = Int
	    required = false
	    default = MaxInt
        "--maxtime_lb_rlxlag"
	    help = "maximum time, in seconds, to compute LB by Lagrangean Relaxation"
            arg_type = Int
	    required = false
	    default = MaxInt
	"--maxtime_lb_colgen"
            help = "maximum time, in seconds, to compute LB by Column Generation"
            arg_type = Int
	    required = false
	    default =  MaxInt
        "--maxtime_ub_lp"
            help = "maximum time, in seconds, to compute UB by LP Rounding"
            arg_type = Int
	    required = false
	    default = MaxInt
	"--maxtime_ub_rlxlag"
            help = "maximum time, in seconds, to compute UB by Lagrangean Relaxation"
            arg_type = Int
	    required = false
	    default = MaxInt
	"--maxtime_ub_colgen"
            help = "maximum time, in seconds, to compute UB by Column Generation"
            arg_type = Int
	    required = false
	    default = MaxInt
	"--maxtime_ilp"
            help = "maximum time, in seconds, to compute UB by Branch and Cut"
            arg_type = Int
	    required = false
	    default =  MaxInt
        "--ra"
            help = "academic number"
	    arg_type = Int
	    required = false
	    default =  999999
    end
    return parse_args(s)
end

function getparameters()
    inputfilename=""
    seednumber=1
    maxtime_lb_lp=MaxInt
    maxtime_lb_rlxlag=MaxInt
    maxtime_lb_colgen=MaxInt
    maxtime_ub_lp=MaxInt
    maxtime_ub_rlxlag=MaxInt
    maxtime_ub_colgen=MaxInt
    maxtime_ilp=MaxInt
    ra=999999
	 
    parsed_args = parse_commandline()
   
    for (arg,val) in parsed_args
	if "$arg"=="seednumber"
	   seednumber=parse(Int,"$val")
	end
	if "$arg"=="inputfilename"
	   inputfilename="$val"
	   if (!isfile(inputfilename))
	      error("Could not find file "*inputfilename)
	   end
	end
	if "$arg"=="ra"
	   ra=parse(IntType,"$val")
	   if (ra<0)
	      error("RA - Academic Number must be a positive integer")
	   end
	end
	if "$arg"=="maxtime_lb_lp"
	   maxtime_lb_lp=parse(IntType,"$val")
	   if (maxtime_lb_lp<0)
	      error("maxtime_lb_lp must be a positive integer (time in seconds)")
	   end
	end
	if "$arg"=="maxtime_lb_rlxlag"
	   maxtime_lb_rlxlag=parse(IntType,"$val")
	   if (maxtime_lb_rlxlag<0)
	      error("maxtime_lb_rlxlag must be a positive integer (time in seconds)")
	   end
	end
	if "$arg"=="maxtime_lb_colgen"
	   maxtime_lb_colgen=parse(IntType,"$val")
	   if (maxtime_lb_colgen<0)
	      error("maxtime_lb_colgen must be a positive integer (time in seconds)")
	   end
	end
	if "$arg"=="maxtime_ub_lp"
	   maxtime_ub_lp=parse(IntType,"$val")
	   if (maxtime_ub_lp<0)
	      error("maxtime_ub_lp must be a positive integer (time in seconds)")
	   end
	end
	if "$arg"=="maxtime_ub_rlxlag"
	   maxtime_ub_rlxlag=parse(IntType,"$val")
	   if (maxtime_ub_rlxlag<0)
	      error("maxtime_ub_rlxlag must be a positive integer (time in seconds)")
	   end
	end
	if "$arg"=="maxtime_ub_colgen"
	   maxtime_ub_colgen=parse(IntType,"$val")
	   if (maxtime_ub_colgen<0)
	      error("maxtime_ub_colgen must be a positive integer (time in seconds)")
	   end
	end
	if "$arg"=="maxtime_ilp"
	   maxtime_ilp=parse(IntType,"$val")
	   if (maxtime_ilp<0)
	      error("maxtime_ilp must be a positive integer (time in seconds)")
	   end
	end
    end

    return(inputfilename,	# input file 
	 seednumber,		# Randomized steps become deterministic for fixed seed
	 maxtime_lb_lp,		# Max. time to compute LB via LP (may use cutting planes)
	 maxtime_lb_rlxlag,	# Max. time to compute LB by Lagrangean Relaxation
	 maxtime_lb_colgen,	# Max. time to compute LB via LP column generation
	 maxtime_ub_lp,		# Max. time to compute UB via LP rounding
	 maxtime_ub_rlxlag,	# Max. time to compute UB via Lagrangean Relaxation
	 maxtime_ub_colgen,	# Max. time to compute UB via LP column generation
	 maxtime_ilp,		# Max. time to compute LB and UB via exact Branch and Cut
	 ra)			# Academic Number
end



