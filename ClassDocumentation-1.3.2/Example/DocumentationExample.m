function [A,B] = DocumentationExample(aninput,options)
%summary of this function, method, or property
%
% Concise, but complete description of this function and how to use it.
%
% - Topic: A topic — A subtopic – A sub-subtopic
% - Declaration: [A,B] = DocumentationExample(aninput,options)
% - Parameter aninput: double
% - Parameter a: string
% - Parameter b: double
% - Returns A: a double
% - Returns B: a double
arguments
    aninput double
    options.a string = ""
    options.b double = 5
end

A = aninput + optiona.b;
B = aninput - options.b;
