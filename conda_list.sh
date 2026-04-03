for env in $(conda env list | grep '^[^#]' | awk '{print $1}'); do 
    echo "--- Environment: $env ---"; 
    conda list -p $env; 
done
