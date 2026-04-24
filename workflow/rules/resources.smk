import os

def get_resource(wildcards, input, attempt, rule_key, resource_name, java_mem_ratio=0.85):
    """
    Generic resource calculator for threads, memory, time, and partition.
    """
    # 1. Access the specific rule configuration or fallback to empty dict
    res_cfg = config.get("resources", {}).get(rule_key, {})

    # 2. Calculate Memory (needed for both 'mem_mb' and 'partition' logic)
    # Base memory falls back to simple_job_mem_mb if rule is not in 'resources'
    base_mem = res_cfg.get("base_mem_mb", config.get("simple_job_mem_mb", 4000))
    scaling_factor = res_cfg.get("scaling_factor", 0)
    
    # Calculate input size in MB
    input_size = 0
    if input:
        if hasattr(input, "size_mb"):
            input_size = input.size_mb
        elif isinstance(input, (list, tuple, dict)):
            # Sum sizes if multiple inputs exist
            for item in input:
                input_size += getattr(item, "size_mb", 0)
    
    # Dynamic Memory Calculation
    mem_mb = (base_mem + (input_size * scaling_factor)) * attempt
    
    # Apply Global Cap
    max_cap = config.get("max_mem_mb", mem_mb)
    mem_mb = int(min(mem_mb, max_cap))

    # 3. Return the specific resource requested
    if resource_name == "mem_mb":
        return mem_mb

    if resource_name == "java_mem":
        return int(mem_mb * java_mem_ratio)

    if resource_name == "threads":
        return res_cfg.get("threads", config.get("simple_job_threads", 1))

    if resource_name == "time_min":
        return res_cfg.get("time_min", config.get("simple_job_time_min", 60))

    if resource_name == "account":
        return res_cfg.get("account", config.get("account", ""))

    if resource_name == "partition":
        ranges = config.get("partition_ranges", {})
        if not ranges:
            return ""
        
        # Find the smallest partition that fits the memory requirement
        # Sorted by limit value
        sorted_ranges = sorted(ranges.items(), key=lambda x: x[1])
        for p_name, p_limit in sorted_ranges:
            if mem_mb <= p_limit:
                return p_name
        # Fallback to the largest partition if it exceeds all thresholds
        return sorted_ranges[-1][0]

    return 0
