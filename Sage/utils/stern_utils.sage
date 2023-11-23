load('list_sorting_utils.sage')
load('LB_utils.sage')


def log2(x):
    return log(x*1.)/log(2.);
#####################################################


########################################################
def PGE(n,r,ell,Fq,H):

#copy_H = copy(H);

    ok = 1;
    i = 0;
    while (i<(r-ell))&(ok):

        #swap rows so that you have a pivot
        j = r-1-i;
        while(H[j,n-1-i]==0)&(j>=0):
            j = j-1;

        #if no valid element has been found, report failure
        if (H[j,n-1-i]==0):
            ok = 0;

        else: #swap rows
            tmp = H[r-1-i,:];
            H[r-1-i,:] = H[j,:];
            H[j,:] = tmp;

            #scale row so that you have the pivot
            scale_coeff = H[r-1-i,n-1-i]^-1;
            H[r-1-i,:] = scale_coeff*H[r-1-i,:];

            #create zeros
            for v in range(r):
                if v!= i:
                    scale_coeff = H[r-1-v,n-1-i];
                    H[r-1-v,:] = H[r-1-v,:] - scale_coeff*H[r-1-i,:];

#            print(H)
#            print("--------------");
        i+=1;

    return ok, H;

###################################################################


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

###Stern's ISD

def stern_isd(H, Fq, Fq_star, n, k, p, ell, w):
    ok = 0;
    while ok == 0:
        P = Permutations(n).random_element().to_matrix();
        new_H = H*P;

        ok, H_prime = PGE(n,n-k,ell,Fq,new_H);

    E = H_prime[0:ell, 0:k+ell];

    #Split the matrix E into submatrices
    E_left = E[:,0:floor((k+ell)/2)];
    E_right = E[:,floor((k+ell)/2):k+ell];

    #create lists for Stern
    list_all_vectors_left = enumerate_Hamming_sphere(Fq, Fq_star, floor((k+ell)/2), p);
    list_subsums_left = create_list_partial_sums(E_left, list_all_vectors_left, matrix(Fq,1,ell),1);

    list_all_vectors_right = enumerate_Hamming_sphere(Fq, Fq_star, ceil((k+ell)/2), p);
    list_subsums_right = create_list_partial_sums(E_right, list_all_vectors_right, matrix(Fq,1,ell),-1);

    #Merge lists
    partial_sums, solutions_small_instance  = merge_lists(list_subsums_left, list_subsums_right, list_all_vectors_left, list_all_vectors_right, matrix(Fq,1,ell), E, 1);


    #list with found codewords
    X = [];
    for e1 in solutions_small_instance:

        e2 = -e1*H_prime[ell:n-k, 0:k+ell].transpose();
        wt_e2 = n-k-ell-e2.list().count(0);
        if wt_e2 == (w - 2*p):
            found_cw = block_matrix(Fq,1,2,[e1,e2])*P.change_ring(Fq)^-1;
            X.append(found_cw);

    return X;
