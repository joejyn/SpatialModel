function [D,t] = Spatial1D_localmigration_GaussianApprox(N,twoNmu,twoNs,twoNm,T,ND)

%Output 
%D is a cell array that contains the frequency array in each deme -
%it has been constucted so that the rows of this array correspond to
%independent origins/haplotypes, with common ordering across all demes, and
%so the number of rows should be the same for all demes

%Input:
% Here N is the population size in each deme
% ND is the number of demes
% T is the number of generations in the simulations
% twoNmu, twoNs & twoNm are the population scaled  - but scaled by the
% total population.
% ND = 3; % number of demes
% N = 10; % population size
% twoNmu = 10;
% s = 0.05; % selection coefficient

%% Population Parameters
mu = twoNmu/2/N/ND; % mutation coefficient
s = twoNs/2/N/ND;
mig = twoNm/2/N/ND; % migration coefficient

%% Preallocating dataframes
migCell = cell(1,ND); % stores how many migrants per deme
migLRA = cell(1,ND); % stores left right array for respective demes. Positive binornd means go to the left
mutantCount = zeros(ND,T); % counts number of de novo mutants generated per deme for each time
migrantCount = zeros(ND,T); % counts number of migrants generated per deme for each time

%% Creating the initial frequency array, with each array representing one deme
D = cell(1,ND); %deme cell array - each with their own frequency array
x = zeros(1,T);
x(end,1) = 1; %Initial frequency array of the wt (no mutants present yet)

for i = 1:ND % populate deme cell array with initial wt-frequency of 1 @ T = 0
    D{1,i} = x;
end

%% Filling up the frequency array
for t = 2:T %No. of generations
    
    for i = 1:ND %generate mutants and drift selection based on t-1, then add in generated mutants for t
        
        x = D{1,i}; % accesses frequency array in selected deme in for loop
                
        %% Determining the mutant and wt frequencies after selection
        xwt = x(end,t-1); %last row always wild type frequency
        X = x(1:end-1,t-1); % extracts the mutant frequencies (no wild type) from previous time point - hence t-1
        xx = X + xwt*s*X; % new mutant frequencies after selection coefficient is considered 
        xx(xx<0)=0; % removes negative values (set as 0)
        xxK = 1-sum(xx); %Frequency of wild type (wt is the Kth allele)
        if xxK < 0 % in case wild type frequency is negative
            xxK = 0; % set wt frequency to 0
            xx = xx/sum(xx); %resets all mutant frequencies to sum to 1
        end     
        
        %% Sample N number of individuals based on drift, using Gaussian approximation of multinomial sampling
        z = GaussApprox_multinomial(N,xx);
        x(1:end,t) = N*z; % converted frequency to individuals to make migrant and mutation calculations easier
        
        %% Determining the number of de novo mutants
        m = poissrnd(N*mu*xwt); %number of de novo mutants generated
        mutantCount(i,t) = m; %stores number of mutants generated into the array

        if m ~= 0 % when a mutant is generated
            
            xN = zeros(m,T); % create an array row (mutant) and column (time)
            xN(1:end,t) = 1; % 1 individual arising per mutant seeding into above array
            x = cat(1,xN,x); % stacks new array ontop of current population of wild type and existing mutants
            
        end

        %% Adding blank haplotypes to all other demes corresponding to the new de novo mutants generated in this deme
        % This ensures that ordering of the haplotypes in D is the same for each deme
        D{1,i} = x; % assigns generated "frequency" array back into deme! still in population numbers hence ""
        
        for f = 1:ND % fill in rows zeroes to accommodate later migrants in other demes
        
            if f == i % skips deme being run in the bigger for-loop
                continue % skips just this iteration of the for-loop
            end
            
            xOD = D{f}; % accesses "frequency" arrays of other demes
            mutantZ = zeros(m,T); % rows of zero for migrants. Size m becuase m de novo mutants were generated in this deme
            xOD = cat(1,mutantZ,xOD); % stacking new rows on top of other deme x arrays
            D{f} = xOD; % stores new "frequency" arrays into respective demes
                     
        end
        
    end
    
    %% Determining the number of migrants leaving the demes for the current time
    for i = 1:ND

        xxM = D{1,i}; % extracting current x array of current deme iteration
        
        MIG = mig*xxM(:,t-1); % generate expected number of migrants per generation based off migrant coefficient 

        migCell{1,i} = poissrnd(MIG)'; % generate actual number of migrants drawn from poisson distribution based off previous time, t-1
        migrantCount(i,t) = sum(migCell{1,i}); % stores number of migrants generated per deme per generation into array
        
    end
    
    %% Determining if migrants move left or right
    for i = 1:ND
        
        [~,c] = size(migCell{1,i}); % registers size of migrant array to preallocate size of migLR array
        migLR = zeros(2,c); % stores how many going left/ right (row 1 left, row 2 right)
        
        if i == 1
            
            migLR(2,:) = migCell{1,i}; % goes right, since only right should be allowed
            migLRA{1,i} = migLR; % stores into left-right array to be used later
            
        elseif i == ND
            
            migLR(1,:) = migCell{1,i}; % goes left, since only left should be allowed
            migLRA{1,i} = migLR; % stores into left-right array to be used later
            
        else % for the middle demes
            
            migLR(1,:) = binornd(migCell{1,i},0.5); % positive goes left therefore row 1
            migLR(2,:) = migCell{1,i} - migLR(1,:); % deducts left migrants from total migrants (as it can only go right if its not going left)
            migLRA{1,i} = migLR; % stores into left-right array to be used later
            
        end
        
    end
    
    %% Accounting for migrants leaving in the "frequency" array, D
    for i = 1:ND % migrants leaving
        
        migLA = migCell{1,i}; % leaving array! not left-right array (migLRA)
        xMig = D{1,i}; % accesses current "frequency" array of deme
        nMig = xMig(:,t); % working in population numbers instead of frequency
        nMig = nMig - migLA'; % migrants leaving
        
        for g = 1:numel(nMig) % poisson might generate more than there is, leading to negative frequency
            
            if nMig(g,1) < 0 % looks for negative values and sets to zero 
                nMig(g,1) = 0;       
            end
            
        end
        
        xMig(:,t) = nMig; % storing population after migrants leave in transition deme array !COLUMN VECTOR
        D{1,i} = xMig; % stores "frequency" array back in to deme cell array
        
    end
    
    %% Accounting for migrants entering in the "frequency" array, D
    for i = 1:ND 
        
        if i == 1 % for first deme
            
            nB = D{1,i}; % n before accounting for entering migrants
            nEA = migLRA{1,i+1}; % n entering array !ROW VECTOR
            nE = nEA(1,:)'; % n entering (only left migrations therefore row 1)!COLUMN VECTOR
            
            nB(:,t) = nB(:,t) + nE; % adds entering array into existing array
            D{1,i} = nB; % stores back into deme cell array
              
        elseif i == ND % for last deme
            
            nB = D{1,i}; % n before
            nEA = migLRA{1,i-1}; % n entering array !ROW VECTOR
            nE = nEA(2,:)'; % n entering (only right migrations therefore row 2)!COLUMN VECTOR
            
            nB(:,t) = nB(:,t) + nE; % adds entering array into existing array
            D{1,i} = nB; % stores back into deme cell array
            
        else % for the middle demes
            
            nB = D{1,i}; % n before !COLUMN VECTOR
            nEAR = migLRA{1,i+1}; % n entering array from right deme !ROW VECTOR
            nER = nEAR(1,:)'; % n entering (only left migration therefore row 1)!COLUMN VECTOR
            nEAL = migLRA{1,i-1}; % n entering array from left deme !ROW VECTOR
            nEL = nEAL(2,:)'; % n entering (only right migration therefore row 2)!COLUMN VECTOR
            
            nB(:,t) = nB(:,t) + nER + nEL; % adds entering array into existing array
            D{1,i} = nB; % stores back into deme cell array
            
        end
        
    end
    
    %% Setting all frequencies to sum up to 1 after migration (convert from population into frequency)
    for i = 1:ND
    
        xF = D{1,i}; % accesses deme cell array
        xF(:,t) = xF(:,t)/sum(xF(:,t)); % population divided by total gives frequency
        D{1,i} = xF; % stores back into cell array
    
    end
    
end

%% Get a vector of generations to plot the simulation graph
t = 1:T;