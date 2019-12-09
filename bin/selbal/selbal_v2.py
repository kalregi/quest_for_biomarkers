from import_data import *
from data_files import *

def balance(X_i, X_j):
    if len(X_i) != len(X_j):
        raise("Columns' length not equal")
    k = len(X_i)
    B = math.sqrt(1/2)*(np.log(X_i) - np.log(X_j))
    return B

def get_AUC(x_train, y_train):

    y_score = sigmoid(x_train)

    #ROC görbe számolása
    fpr, tpr, _ = roc_curve(y_train, y_score)
    roc_auc = auc(fpr, tpr)
    return roc_auc


def opt_bal_for_two(train_table, train_target, test_table, test_target):
    Y_train = [1 if v == "CRC" else 0 for v in train_target]
    Y_test = [1 if v == "CRC" else 0 for v in test_target]
    
    max_col1 = train_table.columns[0]
    max_col2 = train_table.columns[1]
    
    max_auc = 0
    for col1 in train_table.columns:
        for col2 in train_table.columns:
            X_train = np.array(train_table[[col1,col2]])
            B_train = balance(X_train[:,0], X_train[:,1])
            
            X_test = np.array(test_table[[col1,col2]])
            B_test = balance(X_test[:,0], X_test[:,1])
            
            roc_auc = get_AUC(B_train.reshape(-1, 1), Y_train, B_test.reshape(-1, 1), Y_test)
            if roc_auc > max_auc:
                max_auc = roc_auc
                max_col1 = col1
                max_col2 = col2
    return max_col1, max_col2, max_auc

def balance_for_sets(pos_set, neg_set):
    pos_k = pos_set.shape[1]
    neg_k = neg_set.shape[1]
    
    pos_bal = np.zeros((pos_set.shape[0]))
    for k in range(pos_k):
        pos_bal = np.add(pos_bal, np.log(pos_set[:,k]))
    pos_bal = pos_bal * 1/(pos_k)
    
    neg_bal = np.zeros((neg_set.shape[0]))
    for k in range(neg_k):
        neg_bal = neg_bal+np.log(neg_set[:,k])
    neg_bal = neg_bal * 1/(neg_k)
    
    B = pos_bal - neg_bal
    return B

def test_assoc(train_table, train_target, test_table, test_target, pos_set_col_names, neg_set_col_names):
    Y_train = [1 if v == "CRC" else 0 for v in train_target]
    
    pos_set_train = np.array(train_table[pos_set_col_names])
    neg_set_train = np.array(train_table[neg_set_col_names])

    B_train = balance_for_sets(pos_set_train, neg_set_train)  
    
    Y_test = [1 if v == "CRC" else 0 for v in test_target]
    
    pos_set_test = np.array(test_table[pos_set_col_names])
    neg_set_test = np.array(test_table[neg_set_col_names])
    B_test = balance_for_sets(pos_set_test, neg_set_test)

    roc_auc = get_AUC(B_train.reshape(-1, 1),Y_train, B_test.reshape(-1, 1),Y_test)
    return roc_auc

def adding_new_comp(train_table, train_target, test_table, test_target, pos_set_names, neg_set_names, max_auc):
    
    columns = [col for col in train_table.columns if col not in pos_set_names 
               and col not in neg_set_names]
    max_auc_pos = 0
    max_auc_neg = 0
    which = None

    for col in columns:
        
        tmp_pos_set = []
        tmp_pos_set.extend(pos_set_names)
        tmp_pos_set.append(col)
        pos_auc = test_assoc(train_table, train_target, test_table, test_target, tmp_pos_set, neg_set_names)
        
        tmp_neg_set = []
        tmp_neg_set.extend(neg_set_names)
        tmp_neg_set.append(col)
        neg_auc = test_assoc(train_table, train_target, test_table, test_target, pos_set_names, tmp_neg_set)

        if pos_auc > neg_auc:
            if pos_auc > max_auc_pos:
                max_auc_pos = pos_auc
                pos_col = col
                which = "POS"
        else:
            if neg_auc > max_auc_neg:
                max_auc_neg = neg_auc
                neg_col = col
                which = "NEG"
    
    if which == "POS":
        pos_set_names.append(pos_col)
        max_auc_2 = max_auc_pos
    elif which == "NEG":
        neg_set_names.append(neg_col)
        max_auc_2 = max_auc_neg
        
    return pos_set_names, neg_set_names, max_auc_2
    
def opt_bal(train_table, train_target, test_table, test_target, threshold = 0, max_comp = 20):
    print("Optimal balance for 2...")
    max_col1, max_col2, max_auc = opt_bal_for_two(train_table, train_target, test_table, test_target)
    
    pos_set_names = [max_col1]
    neg_set_names = [max_col2]
    
    auc_list = [max_auc]
    
    comp_size = len(pos_set_names) + len(neg_set_names)
    
    crit = True
    print("Selecting further features...")
    while(crit):
        #print("Choosing feature", comp_size+1)
        tmp_pos, tmp_neg, tmp_auc = adding_new_comp(train_table, train_target, test_table, test_target,
                                                    pos_set_names, neg_set_names, max_auc)
        auc_groth = tmp_auc-max_auc
        #print("new auc", tmp_auc)
        if auc_groth > threshold:
            pos_set_names = tmp_pos
            neg_set_names = tmp_neg
            auc_list.append(tmp_auc)
            max_auc = tmp_auc
            comp_size = comp_size + 1
            crit = max_comp > comp_size
        else:
            crit = False
    return pos_set_names, neg_set_names, max_auc, auc_list, comp_size
def test_model(x_train, x_test, y_train, y_test, pos_set_names, neg_set_names):
    
    Y = [1 if v == "CRC" else 0 for v in y_train]
    pos_set = np.array(x_train[pos_set_names])
    neg_set = np.array(x_train[neg_set_names])
    
    B = balance_for_sets(pos_set, neg_set)

    x_train_2 = B.reshape(-1, 1)
    y_train_2 = Y
    
    logisticRegr = LogisticRegression(solver= 'liblinear')
    logisticRegr.fit(x_train_2, y_train_2)
    
    Y_test = [1 if v == "CRC" else 0 for v in y_test]
    pos_set_test = np.array(x_test[pos_set_names])
    neg_set_test = np.array(x_test[neg_set_names])
    
    B_test = balance_for_sets(pos_set_test, neg_set_test)
    
    x_test_2 = B_test.reshape(-1, 1)
    y_test_2 = Y_test

    #Konfidencia scoreok a teszt halmazhoz
    y_score = logisticRegr.decision_function(x_test_2)

    #ROC görbe számolása
    fpr, tpr, _ = roc_curve(y_test_2, y_score)
    roc_auc = auc(fpr, tpr)
    
    plt.figure()
    lw = 2
    plt.plot(fpr, tpr, color='darkorange',
             lw=lw, label='ROC curve (area = %0.2f)' % roc_auc)
    plt.plot([0, 1], [0, 1], color='navy', lw=lw, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title( 'ROC curve')
    plt.legend(loc="lower right")

    return roc_auc

def plot_cv_auc(cross_validation_result):
    auc_lists = []
    pos_set_names = {}
    neg_set_names = {}
    max_comp = 0
    for result in cross_validation_result:
        auc_lists.append(result["auc_list"])
        for name in result["pos_set_names"]:
            if name in pos_set_names:
                pos_set_names[name] = pos_set_names[name] + 1
            else:
                pos_set_names[name] = 1
        for name in result["neg_set_names"]:
            if name in neg_set_names:
                neg_set_names[name] = neg_set_names[name] + 1
            else:
                neg_set_names[name] = 1
        if max_comp < result["comp_size"]:
            max_comp = result["comp_size"]
    
    auc_df = pd.DataFrame(auc_lists, columns = range(2,max_comp+1))
    plt.figure()
    sns.boxplot(data=auc_df)
    return auc_df

def map_selected_cols(cross_validation_result):
    pos_df_list = []
    neg_df_list = []
    fold_num = 0
    for result in cross_validation_result:
        order = 2
        for name in result["pos_set_names"][0:5]:
            pos_df_list.append([fold_num, order, name])
            order = order + 1
        order = 2
        for name in result["neg_set_names"][0:5]:
            neg_df_list.append([fold_num, order, name])
            order = order + 1
        fold_num = fold_num + 1
    pos_df = pd.DataFrame(pos_df_list, columns=["Fold", "Order", "Name"])
    neg_df = pd.DataFrame(neg_df_list, columns=["Fold", "Order", "Name"])
    
    pivot_count_df_neg = pd.DataFrame(neg_df.pivot(index="Name", columns="Fold", values="Order"))
    pivot_count_df_poz = pd.DataFrame(pos_df.pivot(index="Name", columns="Fold", values="Order"))
    joined = pivot_count_df_poz.join(pivot_count_df_neg, lsuffix="_POZ", rsuffix="_NEG")
    
    def poz_neg(row):
        for i in range(10):
            try:
                poz_col_name = str(i) + "_POZ"
                neg_col_name = str(i) + "_NEG"
                if row[poz_col_name] == row[poz_col_name]:
                    if row[neg_col_name] != row[neg_col_name]:
                        row["Fold_"+str(i)] = 1 #POZ
                    else:
                        print("Somethong went wrong, a column can not be in both sets!")
                        raise
                else:
                    if row[neg_col_name] == row[neg_col_name]:
                            row["Fold_"+str(i)] = -1 #NEG
                    else:
                        row["Fold_"+str(i)] = 0 #None
            except:
                pass
        return row

    joined_2 = pd.DataFrame(joined)
    joined_2 = joined_2.apply(poz_neg, axis=1)
    
    keep = []
    for col in joined_2.columns:
        if "Fold_" in col:
            keep.append(col)
    joined_3 = joined_2[keep]
    
    cols = joined_3.columns
    joined_3['Occurance'] = joined_3[cols].gt(0).sum(axis=1) + joined_3[cols].lt(0).sum(axis=1)
    #joined_3["Occurance"] = joined_3.T.count_nonzero()
    joined_3 = joined_3.sort_values("Occurance", ascending=False)
    for_plot = joined_3.drop(columns=["Occurance"])
    
    plt.figure()
    sns.heatmap(for_plot, cmap = sns.color_palette("coolwarm", 7))
    #sns.palplot(sns.color_palette("coolwarm", 7))
    #cmap = sns.diverging_palette(220, 10, as_cmap=True)
    
    return joined_3

def cv_with_test(dataset_name, method, rank, folds = 10, threshold = 0, max_comp = 20, percent=30, table_repl=None, grouping_ser=None):
    start = time.time()

    # Find path to dataset file
    if method == "Kraken":
        table_files = files_kraken
        rank = rank.upper()
    else:
        if method == "Metaphlan":
            table_files = files_metaphlan
            rank = rank.lower()
        else:
            print("Wrong method! Method should be 'Kraken' or 'Metaphlan'")
            raise
    
    
    # Get data
    if table_repl is None and grouping_ser is None:
    	table_repl, grouping_ser = get_table_and_grouping(dataset_name, table_files, rank = rank, percent=percent)
    
    # Split to train and test data
    X_train, X_test, y_train, y_test = train_test_split(table_repl, 
     grouping_ser, test_size=0.2)
    
    # Model selection on train data, with validation
    kf = KFold(n_splits=folds)
    i = 1
    
    cross_validation_result = []

    for train_index, test_index in kf.split(X_train):
        print("Number of fold:", i)
        X_train_k, X_test_k = X_train.iloc[train_index], X_train.iloc[test_index]
        y_train_k, y_test_k = y_train.iloc[train_index], y_train.iloc[test_index]

        pos_set_names, neg_set_names, max_auc, auc_list, comp_size = opt_bal(X_train_k, y_train_k, X_test_k, y_test_k,
                                                              threshold=-1, max_comp=max_comp)
        result = {}
        result["pos_set_names"] = pos_set_names
        result["neg_set_names"] = neg_set_names
        result["max_auc"] = max_auc
        result["auc_list"] = auc_list
        result["comp_size"] = comp_size
        cross_validation_result.append(result)
        i = i+1

    end = time.time()
    print("Elapsed time while CV:", end - start)
    
    #Plots for model selection
    auc_df = plot_cv_auc(cross_validation_result)
    heamap_df = map_selected_cols(cross_validation_result)
    
    # Select comp size (max mean)
    a_l = list(auc_df.mean())
    max_auc_comp_size = a_l.index(max(a_l))+2
    print("Selected cmop size (max mean):", max_auc_comp_size)
    
    # Select comp size (define the optimal number of variables included in the
    #balance as the lowest number whose mean AUC is within 1 standard error of
    # the maximum mean AUC)
    
    # Standard error of mean 
    #np.std(auc_df)/np.sqrt(19)
    # BUT WHERE IS MAXIUMUM?!?!
        # Well I dont care it won't be good anyway....
    
    
    # Train on selected comp size
    X_train_mod, X_test_mod, y_train_mod, y_test_mod = train_test_split(X_train, 
     y_train, test_size=0.2)
    
    pos_set_names, neg_set_names, max_auc, auc_list, comp_size = opt_bal(X_train_mod, 
        y_train_mod, X_test_mod, y_test_mod, threshold, max_auc_comp_size)
    
    print("Positive set:", pos_set_names)
    print("Negative set:", neg_set_names)
    # Training on training data with selected model, than test it on test data
    test_model(X_train, X_test, y_train, y_test, pos_set_names, neg_set_names)
    end2 = time.time()
    print("Elapsed time whole:", end2 - start)
        