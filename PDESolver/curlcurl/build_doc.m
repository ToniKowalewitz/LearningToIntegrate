try
  doc.build;
catch ME
  fprintf('Documentation build failed with message: \n%s\n', ME.message);
  exit(1);
end

exit(0);
