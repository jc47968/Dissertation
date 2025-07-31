# Test script to find the string multiplication error

# Test 1: This should work
message(paste(rep("=", 80), collapse=""))

# Test 2: This should fail 
# message("=" * 80)

# Test 3: Check if there are any hidden characters
test_string <- "="
print(class(test_string))
print(length(test_string))

# Test 4: Manual check
message("Testing message with equals")
message("================================================================================")

message("Debug test completed successfully")