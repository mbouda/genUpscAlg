function testSet=formTestSet()

%% single-layer length laterals
    %one lateral
    testSet(1).parents=[0 1 2 3 3 4 5 6 8]';
    testSet(1).inLayer=[1 2 3 3 3 4 4 5 6]';
    
    %two laterals
    testSet(2).parents=[0 1 2 3 3 4 5 5 6 7 8 9]';
    testSet(2).inLayer=[1 2 3 3 3 4 3 3 5 4 4 6]';

    %three laterals
    testSet(3).parents=[0 1 2 3 3 4 5 5 6 7 7 8 9 10 11]';
    testSet(3).inLayer=[1 2 3 3 3 4 3 3 5 3 3 4 6 4 4]';

    %four laterals
    testSet(4).parents=[0 1 2 3 3 4 5 5 6 7 7 8 9 10 10 11 14 15]';
    testSet(4).inLayer=[1 2 3 3 3 4 3 3 5 3 3 4 6 3 3 4 4 4]';


%% multiple-length laterals

    %three 1-layer, three 2-layer
    testSet(5).parents=[0 1 2 3 3 4 5 5 6 7 7 8 9 10 10 11 12 13 14 14 15 16 19 19 20 21 23 24]';
    testSet(5).inLayer=[1 2 3 3 3 4 3 3 5 3 3 4 6 3 3 4 5 7 3 3 4 5 3 3 4 5 4 4]';

%% unclosed junction with targ basipetal or acropetal to one leg:
    %subs Int going up:
    testSet(6).parents=[0 1 2 3 3 4 5 6 7]';
    testSet(6).inLayer=[1 2 3 3 3 4 4 5 5]';
    
    %subs Int going down:
    testSet(7).parents=[0 1 2 3 3 4 5 6 7 8 9]';
    testSet(7).inLayer=[1 2 3 3 3 4 4 5 5 6 6]';

    %at connUD
    testSet(8).parents=[0 1 2 3 3 4 5 6 7 8 9]';
    testSet(8).inLayer=[1 2 3 4 4 5 5 6 6 7 7]';
    
    %at collar, in same layer
    testSet(10).parents=[0 1 1 2 3 4 5 6 7 8 9]';
    testSet(10).inLayer=[1 1 1 2 2 3 3 4 4 5 5]';
    
    %at collar, across to next layer
    testSet(11).parents=[0 1 1 2 3 4 5 6 7 8 9]';
    testSet(11).inLayer=[1 2 2 3 3 4 4 5 5 6 6]';

%% many long laterals
    %three 5-layer ones
    testSet(9).parents=[0 1 2 3 3 4 5 5 6 7 7 8 9 10 11 12 13 14 15 16 17 18 19 20 22 23 24 25 26]';
    testSet(9).inLayer=[1 2 3 3 3 4 3 3 5 3 3 4 6 4 4 5 7 5 5 6 8 6 6 7 7 7 8 8 8]';

%% unclosed junction with targ acropetal to both legs:

%     %subs Int going  down:
%     testSet(10).parents=[0 1 2 3 3 4 5 6 7 8 9 10 11]';
%     testSet(10).inLayer=[1 2 3 3 3 4 4 5 5 6 6 7 7]';

%at connUD: maybe can't happen with only 2 legs? 
        %because approach is to bring psiL(D+1) up one leg & solve for 2 d.o.f. in other.
        
%% Add more to broaden generality of tests
end