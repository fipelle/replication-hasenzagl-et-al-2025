cd models;
rm ./baseline_oos/data/US_local.xlsx; cp ./baseline_iis/data/US_local.xlsx ./baseline_oos/data/;
rm ./restricted_iis/data/US_local.xlsx; cp ./baseline_iis/data/US_local.xlsx ./restricted_iis/data/;
rm ./restricted_oos/data/US_local.xlsx; cp ./baseline_iis/data/US_local.xlsx ./restricted_oos/data/;
