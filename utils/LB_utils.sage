

def log2(x):
    return log(x*1.)/log(2.);
#####################################################à    
##generate restricted vector with desired weight
def gen_restricted_vector(Fq,E,n,w):
    
    e = vector(Fq,n);
    P = Permutations(range(n)).random_element();
    for i in range(w):
        value = E.random_element();
        e[P[i]] = value;
        
    return e;

##############################


#Create all dispositions of elements from Fq_star, with length num_elements
def all_dispositions(Fq,Fq_star,num_elements):
    
    q = Fq.cardinality();

    all_values = []; #list with all dispositions

    Z = IntegerRing();

    #to generate the elements, we consider all i from 0 to (q-1)^num_elements, and convert each i into
    #its representation in base (q-1).
    #the coefficients of the representation are used to pick the elements from Fq_star
    ##Example:
    #i = 0 ---> repr = (0,0,....,0) ---> e = (Fq_star[0],...,Fq_star[0])
    #i = 1 ---> repr = (0,0,....,1) ---> e = (Fq_star[0],...,Fq_star[1])

    for i in range((q-1)^num_elements): #iterate from 0 to

        repr_i_base_q_minus_1 = Z(str(i)).str(base=(q-1)); #coefficients over base (q-1); note that a has variable size
        max_pos = len(repr_i_base_q_minus_1); #length of i-th representation; we place these coefficients on the rightside

        e = vector(Fq,num_elements);  #vector that contains the disposition associated to i

        #Set the last max_pos entries of e
        for j in range(max_pos):
            e[num_elements-max_pos+j] = Fq_star[int(repr_i_base_q_minus_1[j])];

        #Set the first num_elements - max_pos entries of a to Fq_star[0]
        for u in range(num_elements-max_pos):
            e[u] = Fq_star[0];
            
        #Append created vector and start with a new value of i
        all_values.append(e);
        
    return all_values;     

#######################################################

def create_list_partial_sums(E_sub, list_vectors, target_s, coeff_val):

    list_sums = [];
    for e in list_vectors:
        s = target_s + coeff_val*e*E_sub.transpose();
        list_sums.append(s);

    return list_sums;
#######################################################

#generate all vectors with Hamming weight small_w, length small_n, defined over Fq
#Fq_star must be the multiplicative group of Fq

def enumerate_Hamming_sphere(Fq,Fq_star,small_n,small_w):


    list_of_vectors =  [];#output

    #Generate all possible supports
    all_supports = Combinations(range(small_n),small_w);
    

    if len(Fq_star)>1: #execute if q >= 3

        #Generate all possible dispositions for non null coefficients
        all_values = all_dispositions(Fq,Fq_star,small_w);

        for pos in all_supports: #go through all supports; to each support, we consider all dispositions
            for values in all_values:

                #creating vector with support equal to pos and elements from values
                e = vector(Fq,small_n);

                for u in range(small_w):
                    e[pos[u]] = values[u];

                list_of_vectors.append(matrix(e));

    else: #only if q = 2

        for pos in all_supports:

            e = vector(Fq,small_n);
            
            for u in range(small_w):
                e[pos[u]] = Fq_star[0];
    
            list_of_vectors.append(matrix(e));
    
    return list_of_vectors;

#################################################à        
##Remember to pick n and w as multiples of 4

def prepare_lists(Fq,H,w,E_prime,initial_pos, final_pos, L,b, list_target_syndromes):
    
    initial_lists_of_syndromes = [];
    initial_lists_of_vectors = [];
    
    n = H.ncols();
    small_n = round(n/(2^b));
    small_w = round(w/(2^b));
    
    overall_size = binomial(small_n,small_w)*len(E_prime)^small_w;
    
    #If L is too large, generate all elements and then select only some of them
    if L > 0.01*overall_size:
        
        knapsack_sums = all_subknapsacks(Fq,E_prime,small_n,small_w);

        for i in range(2^b):

            this_H = H[initial_pos: final_pos ,i*small_n:(i+1)*small_n];

            this_list_of_syndromes = [];
            this_list_of_vectors = [];
            
            P = Permutations(range(overall_size)).random_element();
            
            if (i%2) == 0:
                for j in range(L):
                 #   print(j);
                    this_vector = knapsack_sums[P[j]];
                 #   print(this_vector);
                    s = this_vector*this_H.transpose();
                    this_list_of_syndromes.append(s);
                    this_list_of_vectors.append(this_vector);
                    
            else:
                
                target_s = list_target_syndromes[floor(i/2)][0,initial_pos:final_pos];
                for j in range(L):
                    this_vector = knapsack_sums[P[j]];
                    s = target_s - this_vector*this_H.transpose();
                    this_list_of_syndromes.append(s);
                    this_list_of_vectors.append(this_vector);

            initial_lists_of_syndromes.append(this_list_of_syndromes);   
            initial_lists_of_vectors.append(this_list_of_vectors);
        
    else:
        
        #Do not generate the whole list, but extract some random elements
        #generate random lists
        
        pos_list = Combinations(range(small_n),small_w);
        for i in range(2^b):
        
            this_H = H[initial_pos:final_pos, i*small_n:(i+1)*small_n];
            this_list_of_syndromes = [];
            this_list_of_vectors = [];
            
         
            
            if (i%2) == 0:
                for j in range(L):                
                    
                    this_vector = matrix(Fq,1,small_n);
                    this_supp = pos_list.random_element();

                    for x in this_supp:
                        val = E_prime.random_element();
                        this_vector[0,x] = val;
                        
                    
                    s = this_vector*this_H.transpose();
                    this_list_of_syndromes.append(s);
                    this_list_of_vectors.append(this_vector);
            else:
                
                target_s = list_target_syndromes[floor(i/2)][0,initial_pos:final_pos];
                for j in range(L):                
                    
                    this_vector = matrix(Fq,1,small_n);
                    this_supp = pos_list.random_element();

                    for x in this_supp:
                        val = E_prime.random_element();
                        this_vector[0,x] = val;
                        
                    
                    s = target_s - this_vector*this_H.transpose();
                    this_list_of_syndromes.append(s);
                    this_list_of_vectors.append(this_vector);
                
            initial_lists_of_syndromes.append(this_list_of_syndromes);   
            initial_lists_of_vectors.append(this_list_of_vectors);
        
          
            
    return initial_lists_of_syndromes, initial_lists_of_vectors;

################################################

#merge lists for Stern's ISD

def merge_lists(list_syndromes_0,list_syndromes_1, list_vectors_0,list_vectors_1, target_s, sub_H, coeff_val):
    
        
    indexes = colliding_indexes(list_syndromes_0, list_syndromes_1);    

    new_list_syndromes = [];
    new_list_vectors = [];
    for i in range(len(indexes)):

        new_vector = block_matrix(Fq,1,2,[list_vectors_0[indexes[i][0]], list_vectors_1[indexes[i][1]]]);

        new_s = target_s + coeff_val*new_vector*sub_H.transpose();

        new_list_syndromes.append(new_s);
        new_list_vectors.append(new_vector);
    
    return new_list_syndromes, new_list_vectors;


################################################################    

def HJ_W(n,r,b,q,w,v,E_values,ell_values,H,s):

    #Start doing everything
    L = floor(2^v);

    Fq = GF(q);
    E_prime = Set([Fq(x) for x in E_values]);

    small_n = round(n/(2^b));
    

    #Apply random permutation
    P = Permutations(range(n)).random_element();
    new_H = matrix(Fq,r,n);
    for i in range(n):
        for j in range(r):
            new_H[j,i] = H[j,P[i]];
            

    ##Choose target subknapsacks        
    list_target_syndromes = [];
    last_s = s;
    
    for i in range(2^(b-1)-1):
        this_s = random_matrix(Fq,1,r);
        list_target_syndromes.append(this_s);
        last_s = last_s - this_s;
        
    list_target_syndromes.append(last_s);
    
    #Prepare initial lists
    initial_pos = 0;
    final_pos = ell_values[1];
    initial_lists_of_syndromes, initial_lists_of_vectors = prepare_lists(Fq, new_H, w, E_prime, initial_pos, final_pos, L, b, list_target_syndromes);


    #Start merging
    
    initial_pos = ell_values[1];
    final_pos = ell_values[1]+ell_values[2];
  #  print(list_target_syndromes);
    
    for tree_level in range(0,b-1):
#        print(tree_level);

        
        final_list_syndromes = [];
        final_list_vectors = [];
        final_list_target_syndromes = [];
        
    #    print("Doing level = "+str(tree_level));
        
    #    print("=== # of targets = "+str(len(list_target_syndromes)));
    #    print("=== Num of input lists = "+str(len(initial_lists_of_syndromes)));    
        
        avg_list_size = 0;
        
        for i in range(2^(b-tree_level-1)):

 #           print("-------",i);
            
            target_H = new_H[initial_pos:final_pos,i*2^(tree_level+1)*small_n:(i+1)*2^(tree_level+1)*small_n];
                            
            if (i%2) == 0:
                coeff_val = 1;
                target_syndrome = matrix(Fq,1,final_pos-initial_pos);
            else:
                coeff_val = -1;
                
                i0 = i;
                i1 = i-1;
              #  print(initial_pos, final_pos);
              #  print(list_target_syndromes[i0][0,initial_pos:final_pos])
               # target_syndrome = list_target_syndromes[i0][0,initial_pos:final_pos] + list_target_syndromes[i1][0,initial_pos:final_pos];
                target_syndrome = list_target_syndromes[i0] + list_target_syndromes[i1];
              #  print(target_syndrome);
                final_list_target_syndromes.append(target_syndrome);
                target_syndrome = target_syndrome[0,initial_pos:final_pos];
                
            
            new_list_syndromes, new_list_vectors = merge_lists(initial_lists_of_syndromes[2*i],initial_lists_of_syndromes[2*i+1], initial_lists_of_vectors[2*i],initial_lists_of_vectors[2*i+1],target_syndrome, target_H.transpose(), coeff_val);
      #      print("done");
            final_list_syndromes.append(new_list_syndromes);
            final_list_vectors.append(new_list_vectors);
            
            
  #          print("done")
        
        
        #emp_size_values.append(avg_list_size);
        
        if tree_level < (b-2):
            initial_pos = final_pos;
            final_pos = final_pos+ell_values[tree_level+3];
 
        #  print(list_target_syndromes);
        initial_lists_of_syndromes = final_list_syndromes;
        initial_lists_of_vectors = final_list_vectors;
        list_target_syndromes = final_list_target_syndromes;
        
    #    print("---------------------------------------");
        
    final_list_syndromes, solutions = merge_lists(final_list_syndromes[0],final_list_syndromes[1], final_list_vectors[0],final_list_vectors[1], matrix(Fq,1,r), new_H.transpose(),1);
    
    
    
    
    #Remove permutation
    
    final_solutions = [];
    Px = [i+1 for i in P];
    Pm = Permutation(Px).to_matrix();
    Pm = Pm.change_ring(Fq);
    Pinv = Pm^-1;
    for e in solutions:
        ne = e*Pinv;
        final_solutions.append(ne);

    
    
    return final_solutions
