function result = loadFile(fileName, environmentName) 
    tmp = load(fileName, environmentName);
    result = tmp.(environmentName);
end