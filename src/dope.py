

def get_dope():
        total = ''

        with open('model.log', 'r') as log:
                for line in log:
                        total+= line

        dopelist = total.split('>> Summary of successfully produced models:')[-1].split('openf___')[0]


        dopeline = dopelist.split('----------------------------------------------------------------------')[-1].split('\n')

         
        results = []
        for model in dopeline:
                results.append(model.split())


        return results
         

def get_multi_dope():

        total = ''

        with open('model.log', 'r') as log:
                for line in log:
                        total+= line

        dopelist = total.split('>> Summary of successfully produced models:')[-1].split('COMPARISON OF SEVERAL 3D STRUCTURES:')[0]


        dopeline = dopelist.split('----------------------------------------------------------------------')[-1].split('\n')

         
        results = []
        for model in dopeline:
                results.append(model.split())


        return results


get_multi_dope()
