import math


# This class defines a particular "quantum architecture",
# meaning the physical features like error rate are defined 
# for each object, and once these are defined, the costs
# of various algorithms are computed for that object
class QuantumArchitecture:
	def __init__(
		self,
		qubits_per_layer=15,
		qubits_per_final_distill=16,
		qubits_distance=12.5,
		distill_error_constant=35,
		distill_gates=10.0,
		distance_error_constant=0.03,
		distill_cycles=2.0*3.0*1.25,
		qubit_error=0.01,
		threshold_error=0.0316,
		state_injection_error=0.05,
		time_mem=0.5,
		max_fact=1000,
		):
		"""Constructor for the architecture of a quantum computer
		
		Args:
			Q_per_distill: An int that gives the number of qubits each distillation requires
			Q_per_final_distill: An int giving the number of qubits the final distillation requires
			Q_distance_constant: Num. of physical qubits = Num. of logical qubits * Q_distance_constant * distance^2
			Distill_gates: The number of gates used in the distillation circuit
			Distill_error_constant: The constant for how much error increases with distillation
			Distill_cycles: The number of cycles for one distillation
			p: The error rate of physical qubits
			pt: The "threshold error", a constant of the code
			pinject: THe state injection error rate
			time_mem: The time memory tradeoff, in [0,1]
			max_fact: The maximum number of state distillation factories
			To-do: Validate these inputs server-side
		"""
		self.Q_per_distill = int(qubits_per_layer)
		self.Q_per_final_distall = int(qubits_per_final_distill)
		self.Q_distance_constant = float(qubits_distance)
		self.Distill_error_constant = int(distill_error_constant)
		self.Distill_gates = float(distill_gates)
		self.Distance_error_constant = float(distance_error_constant)
		self.Distill_cycles = float(distill_cycles)
		self.p = float(qubit_error)
		self.pt = float(threshold_error)
		self.pinject = float(state_injection_error)
		self.time_mem=float(time_mem)
		self.max_fact=int(max_fact)


	def distill(self,PA):
		#Computes the necessary distances of each distillation to achieve
		#a final error rate of PA
		
		#d=array of distances (bottom up)
		circuit_size=self.Distill_cycles*self.Q_per_final_distall
		d=[self.get_distance(circuit_size,self.p/self.pt,PA)]
		pdist=self.pinject
		L=0
		#Determine number of layers
		while pdist>PA:
			pdist=self.Distill_error_constant*(pdist**3)
			L=L+1

		#Compute distances for each layer based on minimum to achieve fidelity in next layer
		for l in range(1,L-1):
			#This accounts for the extra qubits needed with each layer 
			#because each layer needs more qubits, which increases the necessary fidelity
			#This is somewhat approximate: if we take longer, the error rate must be lower
			#This would be by a factor of, say, 2, for one extra cycle? 
			#I assume this is small: possibly a (tedious) fix for v2
			pdist=math.pow(self.Distance_error_constant*math.pow(self.p/self.pt,d[l-1])/float(self.Distance_error_constant),1.0/3.0);
			circuit_size=circuit_size*self.Q_per_distill
			d.append(self.get_distance(circuit_size,self.p/self.pt,pdist))
		return d

	def get_distance(self,circuit_size,p_ratio,pa ):
		#Finds the distance that will have
		#the right error 
		#Since the error rate depends on the distance
		#in a hard-to-invert way, this function
		#iterates through possible distances until it 
		#is high enough
		dold=1;
		d=2;
		while d!=dold:
			dold=d
			pl = pa/(circuit_size*dold)
			d=math.ceil(math.log(pl)/math.log(p_ratio))
		return d

	def distill_time_space(self,d,rate,correct='true'):
		#Outputs t,q,prod
		#Meaning: It takes each layer time t to produce prod states, and requires
		#q logical qubits
		#Idea: Some input rate (i.e., we have that many logical qubits coming in)
		#This computes the average number distilled for the next layer
		#If it's more than Q_per_distill then the time is just the time to perform distillation
		#And then the next state gets at least one state
		#If rate is too low to instantly distill, then we have to wait, so the production becomes 1
		
		#If we have to save qubits then we can't reuse stages that produce a full set
		#So "correct" will move to an optimized paralellization
		overlap=[1]*len(d)
		if correct and rate>self.Q_per_distill:
			i=0;
			while rate>self.Q_per_distill:
				rate=float(rate)/float(self.Q_per_distill)
				i=i+1
			rate=math.floor(rate)
			for j in range(0,i):
				rate=rate*self.Q_per_distill
				overlap[j]=0
			overlap[i-1]=1
		d.reverse()
		t=[]
		q=[]
		prod=[]
		time=0
		for dist in d:
			if rate>self.Q_per_distill:
				time=time+self.Distill_gates*dist 
				q.append(self.Q_per_distill*math.ceil(rate/self.Q_per_distill));
				rate=float(rate)/float(self.Q_per_distill)
			else:
				time=max(time,1)*(self.Q_per_distill/rate)+self.Distill_gates*dist
				rate=1
				q.append(self.Q_per_distill)
			prod.append(rate)
			t.append(time)
		d.reverse()
		t.reverse()
		q.reverse()
		overlap.reverse()
		prod.reverse()
		return {'t':t,'q':q,'prod':prod,'overlap':overlap}
		
	def surface_code_cost(self,T,D,Q):
		#what this needs to do:
		#Compute necessary distillation and distances
		#Compute total time and memory for distillation
		#if the time is too long, recompute distances and try again
		
		#Depth*qubit is total chances for error: that must be small
		PA=min(1.0/T,1.0/(D*Q))
		d=[1.0];
		total_depth=D
		fact=1
		qubits_per_magic_state=self.Q_distance_constant*(d[0]**2)
		while self.Distance_error_constant*math.pow(self.p/self.pt,d[0]) > PA:
	 		d=self.distill(PA)
			full_mem=self.Q_per_distill**len(d)
			sizes=self.distill_time_space(d,self.time_mem*full_mem+1)
			time_per_magic_state=sizes['t'][0]
			qubits_per_magic_state=self.Q_distance_constant*(d[0]**2)
			#Compute total qubits, by ignoring ones that can overlap
			for i in range(1,len(d)-1):
				if sizes['overlap'][i]==1:
					qubits_per_magic_state=qubits_per_magic_state+self.Q_distance_constant*(d[i]**2)
			#How many state factories are necessary
			fact=math.ceil(T*time_per_magic_state/total_depth)
			#How many it can actually have
			fact = min(fact,self.max_fact)
			#Hence, the total time
			total_depth = max(total_depth,round((T*time_per_magic_state)/fact))
			#The actual necessary fidelity
			PA=min(PA,1.0/(total_depth*Q))
		total_qubits=Q*self.Q_distance_constant*(d[0]**2)+qubits_per_magic_state*fact

		return {'depth':total_depth,'qubits':total_qubits,'distances':d,'factories':fact}


	def algorithm_cost(self,algorithm,inputs):
		#algorithm should be a string that matches one of these options
		#inputs should be a dictionary that matches the next function
		#This retuurns a dictionary of {'depth','qubits','distances','factories'}
		if algorithm=="Generic":
			results=self.generic(**inputs)
		elif algorithm=="Grover-generic":
			results=self.grover(**inputs)
		elif algorithm=="Grover-AES":
			results=self.grover_aes(**inputs)
		elif algorithm=="Shor-Zalka":
			results=self.shor_via_zalka(**inputs)
		elif algorithm=="Shor-Markov":
			results=self.shor_via_markov_saeedi(**inputs)
		elif algorithm=="ECDLP":
			results=self.shor_ecdlp(**inputs)
		else:
			results={"failure":1}
		return results



	#All of the following functions compute
	#the requires T-gates, depth, and qubits
	#then use the previous functions to get the costs

	def generic(self,t_gates,depth,qubits):
		return self.surface_code_cost(int(t_gates),int(depth),int(qubits))


	def shor_via_markov_saeedi(self,bit_length):
		N=int(bit_length)
		T=3.351*(N**3)-145*(N**2)+243*N-104;
		CNOT=5.309*(N**3)-39*(N**2)+76*N-62;
		return self.surface_code_cost(T,max(T,CNOT),5*N+1)
			
	def shor_via_estera_hastad(self,bit_length):
		N=int(bit_length)
		T=3.351*(N**3)-145*(N**2)+243*N-104;
		CNOT=5.309*(N**3)-39*(N**2)+76*N-62;	
		T=T/2;
		CNOT=CNOT/2;
		return self.surface_code_cost(T,max(T,CNOT),5*N+1)
		
	def shor_via_zalka(self,bit_length):
		N=int(bit_length)
		return self.surface_code_cost(7*52*(N**3),3*500*(N**2),5*N)
		
	def shor_via_fowler(self,bit_length):
		N=int(bit_length)
		return self.surface_code_cost(280*(N**3),120*(N**3),2*N)
		
	def bbm(N):
		delta=math.pow(8.0/3.0,1.0/3.0)
		beta=math.pow(1.0/3.0,1.0/3.0)
		epsilon=3.0*beta/2.0
		Ne=N*math.log(2)
		X=(2.0/delta+delta*epsilon)*math.pow(Ne,2.0/3.0)*math.pow(math.log(Ne),1.0/3.0)
		T=280*((X/math.log(2))**3)
		D=120*((X/math.log(2))**3)
		Q=2*X/math.log(2)
		R=math.exp((epsilon-beta/2.0)*math.pow(Ne,1.0/3.0)*math.pow(math.log(Ne),2.0/3.0))
		return self.surface_code_cost(T*R,D*R,Q)

	def grover(self,oracle_t_gates,oracle_depth,qubits,search_space):
		oT=int(oracle_t_gates)
		oD=int(oracle_depth)
		Q=int(qubits)
		N=int(search_space)
		#Multiplies the number T-gates/deph
		#per oracle call by the square root of the seach space
		T=round(oT*math.sqrt(N)*2)
		D=round(oD*math.sqrt(N)*2)
		return self.surface_code_cost(T,D,Q)
		
	def grover_aes(self,AES_length):
		N=int(AES_length)
		#from grassl et al.
		if (N==128):
			T=1.19*(2**86)
			Q=2953
			D=1.16*(2**81)
		elif (N==192):
			T=1.81*(2**118)
			Q=4449
			D=1.33*(2**113)
		elif (N==256):
			T=1.41*(2**151)
			Q=6681
			D=1.57*(2**145)
		return self.surface_code_cost(T,D,Q)

		
	def shor_ecdlp(self,prime_length):
		P=int(prime_length)
		T={110:(9.44*(10**9)),
			160:(2.97*(10**10)),
			192:(5.30*(10**10)),
			224:(8.43*(10**10)),
			256:(1.26*(10**11)),
			384:(4.52*(10**11)),
			521:(1.14*(10**12))}
		Q={110:1014,
			160:1466,
			192:1754,
			224:2042,
			256:2330,
			384:3484,
			521:4719}
		D={110:(8.66*(10**9)),
			160:(2.73*(10**10)),
			192:(4.86*(10**10)),
			224:(7.73*(10**10)),
			256:(1.16*(10**11)),
			384:(4.15*(10**11)),
			521:(1.05*(10**12)),
			}
		return self.surface_code_cost(T[P],D[P],Q[P])
	

def runTests():
	test=QuantumArchitecture(
			15,
			16,
			12.5,
			35,
			10.0,
			0.03,
			2.0*3.0*1.25,
			0.01,
			0.0316,
			0.05,
			0.5,
			1000
			)
	#various tests
	#generic
	result_generic=test.algorithm_cost("Generic",{"tgates":20,"depth":30,"qubits":50})
	print("Any algorithm:")
	print(result_generic)
	#generic grover
	results_grover_generic=test.algorithm_cost("Grover-generic",{"oracletgates":20,"oracledepth":30,"qubits":50,"searchspace":1000000000})
	print("Any Grover:")
	print(results_grover_generic)
	#grover AES
	results_grover_AES=test.algorithm_cost("Grover-AES",{"AESlength":192})
	print("Grover AES:")
	print(results_grover_AES)
	#Shor via Zalka
	results_shor_zalka=test.algorithm_cost("Shor-Zalka",{"bitlength":2048})
	print("Shor with Zalka")
	print(results_shor_zalka)
	#Shor via Markov and Saeedi
	results_shor_markov=test.algorithm_cost("Shor-Markov",{"bitlength":2048})
	print("Shor with Markov")
	print(results_shor_markov)
	#ECDLP
	results_ECDLP=test.algorithm_cost("ECDLP",{"primelength":384})
	print("ECDLP:")
	print(results_ECDLP)

#runTests()
