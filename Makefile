CONFIG_FILE = config.mk
include $(CONFIG_FILE)

# Encontra todos os arquivos .txt no diretório de instâncias
# A função $(wildcard ...) expande para uma lista de arquivos que correspondem ao padrão
# A função $(notdir ...) remove o caminho do diretório, deixando apenas o nome do arquivo
# A função $(basename ...) remove a extensão do arquivo
INSTANCES := $(basename $(notdir $(wildcard $(INSTANCES_DIR)/*.txt)))

.PHONY: all test clean

# Regra 'all': executa o problema para cada instância
all: $(foreach instance, $(INSTANCES), $(instance).log)

# Regra de pattern matching para gerar cada arquivo de log
# $@ é o nome do target (ex: grf101.log)
# $* é o nome do target sem a extensão (ex: grf101)
%.log:
	@echo "Executando para a instância: $*"
	julia --project=. trigger_arc_tsp.jl \
		--inputfilename $(INSTANCES_DIR)/$*.txt \
		--seednumber $(SEED_NUMBER) \
		--maxtime_lb_lp $(MAXTIME_LB_LP) \
		--maxtime_lb_rlxlag $(MAXTIME_LB_RLXLAG) \
		--maxtime_lb_colgen $(MAXTIME_LB_COLGEN) \
		--maxtime_ub_lp $(MAXTIME_UB_LP) \
		--maxtime_ub_rlxlag $(MAXTIME_UB_RLXLAG) \
		--maxtime_ub_colgen $(MAXTIME_UB_COLGEN) \
		--maxtime_ilp $(MAXTIME_ILP) \
		--ra $(RA_PARAM) \
		--logfile logs/test_$*.log \
		$(VERBOSE_MODE)

test:
	@echo Execução de instância de teste:
	julia --project=. trigger_arc_tsp.jl \
			--inputfilename $(TEST_INSTANCE) \
			--seednumber 1234 \
			--maxtime_lb_lp 1010 \
			--maxtime_lb_rlxlag 1020 \
			--maxtime_lb_colgen 1030 \
			--maxtime_ub_lp 1040 \
			--maxtime_ub_rlxlag 1050 \
			--maxtime_ub_colgen 1060 \
			--maxtime_ilp 60 \
			--ra 291111 \
			--logfile trigger_arc_tsp.log

# Regra 'clean': remove todos os arquivos de log gerados
clean:
	@echo "Removendo arquivos de log..."
	rm -f $(INSTANCES)/logs/*.log