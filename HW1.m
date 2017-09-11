%% Homework 1. Due before class on 9/5/17


%Jiahui Madelaine Li (jl95)
%Fall 2017, BIOC 470

% Fill in the blank space in this section with code that will add 
% the two numbers regardless of variable type. Hint see the matlab
% functions ischar, isnumeric, and str2num. 
function HW1
    %HW1 function call other functions to solve the problems in homework 1
    
    %your code should work no matter which of these lines is uncommented. 
    x = 3; y = 5; % integers
    %x = '3'; y= '5'; %strings
    %x = 3; y = '5'; %mixed
    
    %problem 1
    stringaddition(x,y);
    
    %problem 2
    problem_2_a = randseq(500);
    disp('problem 2, part 1, random sequence');
    disp(problem_2_a);
    
    problem_2_b = find_ORF(500);
    disp('problem 2, part 2, longest ORF');
    disp(problem_2_b);
    
    problem_2_c = prob_50_same_length(500);
    disp('problem 2, part 3, probability of >50b.p.');
    disp(problem_2_c);
    
    figure(1)
    prob_50_diff_length(500);
    
    %problem 3
    problem_3_a = read_data();
    disp('problem 3, part 1, vector of Cp');
    disp(problem_3_a);
    
    problem_3_b = vector_into_array();
    disp('problem 3, part 2, array of qPCR layout');
    disp(problem_3_b);
    
    figure(2)
    normalization();
end
%% Problem 1 - addition with strings
%your code goes here
function stringaddition(x,y)
    if ischar(x)==1
        x = str2num(x);
    end
    if ischar(y)==1
        y = str2num(y);
    end
    sum = x+y;
    
    %output your answer
    disp ('problem 1')
    disp (sum)
end

%% Problem 2 - our first real biology problem. Open reading frames and nested loops.
%part 1: write a piece of code that creates a random DNA sequence of length
% N (i.e. consisting of the letters ATGC) where we will start with N=500 base pairs (b.p.).
% store the output in a variable
% called rand_seq. Hint: the function randi may be useful. 
% Even if you have access to the bioinformatics toolbox, 
% do not use the builtin function randseq for this part. 
function rand_seq = randseq(N)
    %generate a random DNA sequence of length N
    randnum = randi(4,[1,N]);
    for x = 1:N
        if randnum(x) == 1
            rand_seq(x) = 'A';
        elseif randnum(x) == 2
            rand_seq(x) = 'G';
        elseif randnum(x) == 3
            rand_seq(x) = 'C';
        elseif randnum(x) == 4
            rand_seq(x) = 'T';
        end
    end
end
%part 2: open reading frames (ORFs) are pieces of DNA that can be
% transcribed and translated. They start with a start codon (ATG) and end with a
% stop codon (TAA, TGA, or TAG). Write a piece of code that finds the longest ORF 
% in your seqeunce rand_seq. Hint: see the function strfind.
function longest_ORF = find_ORF(N)
    rand_seq = randseq(N);
    
    %find start codons' positions
    start_ind = strfind(rand_seq, 'ATG');
    %find end codons' positions
    end_ind_unsorted = [strfind(rand_seq, 'TAA'), strfind(rand_seq, 'TGA'), strfind(rand_seq, 'TAG')];
    end_ind = sort(end_ind_unsorted);

    ORF = [];
    
    %for every start codon
    for ind = 1:length(start_ind)
        diff_array = end_ind - start_ind(ind);
        %only find frames that are multiple of 3
        isframe = mod(diff_array,3);
        for ii = 1:length(isframe)
            if isframe (ii) ~=0
                diff_array(ii) = NaN;
            end
        end
        %only find frames that exist (start codon in front of end codon)
        viable = find(diff_array>0);
            if isempty(viable)
                ORF(ind) = 0;
                ORF_start(ind) = 0;
                ORF_end(ind) = 0;
            else
                ORF(ind) = diff_array(viable(1));
                ORF_start(ind) = start_ind(ind);
                ORF_end(ind) = end_ind(viable(1));
            end
    end
    %find longest ORF if there are ORFs starting from this start codon
    longest = max(ORF);
    if longest~=0
        idx = find(ORF == longest);
        start_nuc = ORF_start(idx(1));
        end_nuc = ORF_end(idx(1));
        longest_ORF = rand_seq(start_nuc:end_nuc+2);
    else
        longest_ORF = 'none';
    end
end



%part 3: copy your code in parts 1 and 2 but place it inside a loop that
% runs 1000 times. Use this to determine the probability
% that an sequence of length 500 has an ORF of greater than 50 b.p.
function prob = prob_50_same_length(N)
    total = 0;
    %find longest ORF for 1000 times with a sequence length of N
    for n = 1:1000
        longest_ORF = find_ORF(N);
        ORF_length = length(longest_ORF);
        %count the number of longest ORFs that >50bp
        if ORF_length > 50
            total = total + 1;
        end
    end
    prob = total / 1000;
end


%part 4: copy your code from part 3 but put it inside yet another loop,
% this time over the sequence length N. Plot the probability of having an
% ORF > 50 b.p. as a funciton of the sequence length. 
function prob_50_diff_length(sl)
    prob = [];
    %for all sequence lengths from 1 to sl, find the prob of having ORF>50
    for N = 1:sl
        prob(N) = prob_50_same_length(N);
    end
    plot (1:sl, prob);
    title('problem 2, part 4, probability of >50b.p. vs. sequence length');
    xlabel ('sequence length')
    ylabel ('probability of >50b.p.')
end

%part 5: Make sure your results from part 4 are sensible. What features
% must this curve have (hint: what should be the value when N is small or when
% N is very large? how should the curve change in between?) Make sure your
% plot looks like this. 

%This curve should be sloped upward. The probability should increase as N
%increases. When N is small (<50), probability should equal to 0. When N
%starts to exceed 50, probability should start to increase from 0. 
%When N is very large, probability should approach to 1, but never exceed 1. 
%In between, the probability should increase like a sigmoid function. 


%% problem 3 data input/output and simple analysis

%The file qPCRdata.txt is an actual file that comes from a Roche
%LightCycler qPCR machine. The important columns are the Cp which tells
%you the cycle of amplification and the position which tells you the well
%from the 96 well plate. Each column of the plate has a different gene and
%each row has a different condition. Each gene in done in triplicates so
%columns 1-3 are the same gene, columns 4-6 the same, etc.
%so A1-A3 are gene 1 condition 1, B1-B3 gene 1 condition 2, A4-A6 gene 2
%condition 1, B4-B6 gene2 condition 2 etc. 

% part1: write code to read the Cp data from this file into a vector. You can ignore the last two
% rows with positions beginning with G and H as there were no samples here. 
function vector = read_data()
    %only read Cp data from file
    filename = 'qPCRdata.txt';
    fid = fopen (filename, 'r');
    vector = [];
    for i = 1:74
        line = fgetl(fid);
        text = textscan(line,'%s');
        if i > 2
            vector(i-2) = str2num(text{1}{6});
        end
    end
end

% Part 2: transform this vector into an array representing the layout of
% the plate. e.g. a 6 row, 12 column array should that data(1,1) = Cp from
% A1, data(1,2) = Cp from A2, data(2,1) = Cp from B1 etc. 
function data = vector_into_array()
    vector = read_data();
    data = [];
    %number the vector elements by counting 12 into a row, and starting the
    %next row after12 elements
    for i = 1:72
        row = floor(i/12)+1;
        if mod(i,12) == 0
            row = row -1;
        end
        col = mod(i,12);
        if col ~= 0
            data(row,col) = vector(i);
        else
            data(row,12) = vector(i);
        end
    end
end


% Part 3. The 4th gene in columns 10 - 12 is known as a normalization gene.
% That is, it should not change between conditions and it is used to normalize 
% the expression values for the others. For the other three
% genes, compute their normalized expression in all  conditions, normalized to condition 1. 
% In other words, the fold change between these conditions and condition 1. The
% formula for this is 2^[Cp0 - CpX - (CpN0 - CpNX)] where Cp0 is the Cp for
% the gene in the 1st condition, CpX is the value of Cp in condition X and
% CpN0 and CpNX are the same quantitites for the normalization gene.
% Plot this data in an appropriate way. 
function normalization()
    data = vector_into_array();
    gene = [];
    %take the average for the data from same gene, same condition
    for i = 1:6
        for j = 1:4
            gene(i,j) = (data(i,j*3-2)+data(i,j*3-1)+data(i,j*3))/3;
        end
    end
    
    CpN0 = gene(1,4);
    normed = [];
    for i = 1:4
        for j = 1:6
            Cp0 = gene(1,i);
            CpX = gene(j,i);
            CpNX = gene(j,4);
            normed (j,i) = 2^(Cp0 - CpX - (CpN0 - CpNX));
        end
    end
    %plot the normalized Cp values for each gene against the 6 conditions
    cond = 1:6;
    plot (cond,normed(:,1),cond,normed(:,2),cond,normed(:,3),cond,normed(:,4));
    title('normalized Cp values vs. 6 different conditions');
    xlabel('conditions')
    ylabel('normalized Cp values')
    legend('gene 1','gene 2','gene 3','gene 4')
    
end


%% Challenge problems that extend the above (optional)

% 1. Write a solution to Problem 2 part 2 that doesn't use any loops at
% all. Hint: start by using the built in function bsxfun to make a matrix of all distances
% between start and stop codons. 

% 2. Problem 2, part 4. Use Matlab to compute the exact solution to this
% problem and compare your answer to what you got previously by testing
% many sequences. Plot both on the same set of axes. Hint: to get started 
% think about the following:
% A. How many sequences of length N are there?
% B. How many ways of making an ORF of length N_ORF are there?
% C. For each N_ORF how many ways of position this reading frame in a
% sequence of length N are there?

% 3. Problem 3. Assume that the error in each Cp is the standard deviation
% of the three measurements. Add a section to your code that propogates this
% uncertainty to the final results. Add error bars to your plot. (on
% propagation of error, see, for example:
% https://en.wikipedia.org/wiki/Propagation_of_uncertainty


