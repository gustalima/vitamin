
def get_dope():
        total = ''

        with open('loop_ref.log', 'r') as log:
                for line in log:
                        total+= line

        dopelist = total.split('----------------------------------------')[-1]


         
        results = dopelist.split('\n')


        
        return results
         
